
# Spreadsheet Driven Version
# Date             : Jul 16 2025
#
# Author:          : Jennifer Marr (jennifer.marr@gov.bc.ca)
#
# Purpose:         : Find individual tree crowns from a canopy height model (CHM)
#
# Inputs:          : chm raster and output GDB
#
# Outputs:         : canopy_mask — A binary raster where 1 = canopy (pixels >= 3 m tall) and NoData = everything else. Used internally to clip the watershed results, but also useful to verify the canopy extent visually.
#
#                  : treetops - A binary raster marking which pixels were identified as treetop peaks (local maxima at either the small or large search radius). 
#                               Useful for QA — you can overlay this on the CHM to check whether peaks are being detected in the right places, and tune 
#                               treetop_radius_small/treetop_radius_large accordingly.
#
#                  : crown_polygons — The final deliverable. One polygon per individual crown with fields for area, perimeter, max width, min width, and mean width. The "CrownID" field is just a unique identifier that matches the OBJECTID of the polygon.
#
# History:         : (2026-06-05) Created


"""
Segment canopy crowns from a CHM raster and calculate crown widths.

REQUIRES:
- Spatial Analyst extension

INPUTS:
- CHM raster (canopy height model, in meters)
- Output geodatabase
- Study area polygon
Optional to change:
- Minimum canopy height (m): anything below this height is treated as non-canopy
- Treetop search radius (small, cells): local-maximum window for detecting smaller, tighter crowns
- Treetop search radius (large, cells): local-maximum window for broader crowns; increase for larger trees and flatter tops.
- Minimum crown area (m2): polygons smaller than this threshold are assumed to be noise and removed.

OUTPUTS:
- Crown polygon feature class
- Width metrics as AddMessage: MaxWidth_m, MinWidth_m, MeanWidth_m

NOTES:
- This is an approximate crown delineation workflow.
- Best for reasonably distinct crowns.
- Dense closed-canopy forests may require treetop-guided segmentation.
"""

import arcpy
import os
import math
from arcpy.sa import *

arcpy.management.ClearWorkspaceCache() #clean personal arcpro temp folder
arcpy.env.overwriteOutput = True


# -------------------------------------------------------------------
# USER SETTINGS
# -------------------------------------------------------------------
chm_raster = arcpy.GetParameterAsText(0)
out_gdb = arcpy.GetParameterAsText(1)
study_area = arcpy.GetParameterAsText(2)  # polygon defining the study area boundary

scratch_workspace = "memory"

# Output names
canopy_mask = os.path.join(out_gdb, "canopy_mask")
treetop_raster = os.path.join(out_gdb, "treetops")
treetop_seeds = os.path.join(scratch_workspace, "treetop_seeds")
crown_segments = os.path.join(scratch_workspace, "crown_segments")
crown_polygons_multi = os.path.join(scratch_workspace, "crown_polygons_multi")
crown_polygons = os.path.join(out_gdb, "crown_polygons")

# MIGHT WANT TO TUNE THESE FOR YOUR DATA DEPENDING ON TREE SPECIES AND AREA
min_canopy_height = float(arcpy.GetParameterAsText(3))      # meters; minimum CHM height to count as canopy
treetop_radius_small = int(arcpy.GetParameterAsText(4))     # cells; detects small/narrow crowns
treetop_radius_large = int(arcpy.GetParameterAsText(5))     # cells; detects large/broad crowns (increase for very large trees)
min_crown_area = float(arcpy.GetParameterAsText(6))         # square meters; remove tiny polygons after conversion


# -------------------------------------------------------------------
# HELPER FUNCTIONS
# -------------------------------------------------------------------
def add_field_if_missing(fc, field_name, field_type, **kwargs):
    existing = [f.name for f in arcpy.ListFields(fc)]
    if field_name not in existing:
        arcpy.management.AddField(fc, field_name, field_type, **kwargs)


def polygon_width_metrics(geom):
    """
    Estimate crown width metrics from a polygon geometry.

    Returns:
        area_m2, max_width, min_width, mean_width
    """
    if not geom or geom.area == 0:
        return None, None, None, None

    area = geom.area
    extent = geom.extent

    # Collect all polygon vertices
    pts = []
    part_count = geom.partCount
    for i in range(part_count):
        part = geom.getPart(i)
        for p in part:
            if p:
                pts.append((p.X, p.Y))

    if len(pts) < 2:
        return area, None, None, None

    # Maximum width = max pairwise distance between vertices
    max_dist = 0.0
    for i in range(len(pts)):
        x1, y1 = pts[i]
        for j in range(i + 1, len(pts)):
            x2, y2 = pts[j]
            d = math.hypot(x2 - x1, y2 - y1)
            if d > max_dist:
                max_dist = d

    # Minimum width approximation:
    # use area / (pi/4 * major_axis) scaled as ellipse minor axis estimate
    # ellipse area = pi * a * b / 4 where a,b are diameters
    # so b = 4A / (pi * a)
    min_dist = None
    if max_dist > 0:
        min_dist = (4.0 * area) / (math.pi * max_dist)

    mean_width = None
    if max_dist is not None and min_dist is not None:
        mean_width = (max_dist + min_dist) / 2.0

    return area, max_dist, min_dist, mean_width


# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------
def main():
    arcpy.env.workspace = out_gdb
    arcpy.env.scratchWorkspace = scratch_workspace

    # Check out Spatial Analyst
    arcpy.CheckOutExtension("Spatial")

    # Match processing environment to the CHM so all raster outputs align perfectly.
    desc = arcpy.Describe(chm_raster)
    arcpy.env.extent = chm_raster
    arcpy.env.snapRaster = chm_raster
    arcpy.env.cellSize = desc.meanCellWidth

    chm = Raster(chm_raster)

    # 1. Write anything with a cell height < min_canopy_height as NoData, everything else = 1 (canopy)
    arcpy.AddMessage("Creating canopy mask...")
    mask = Con(chm >= min_canopy_height, 1)
    mask.save(canopy_mask)

    # 2. Find treetops 
    # creates a raster where each cell stores the maximum CHM value in a circular neighborhood.  Twice, once with a small radius and once with a large radius. 
    arcpy.AddMessage("Detecting treetops (multi-scale)...")
    focal_max_small = FocalStatistics(chm, NbrCircle(treetop_radius_small, "CELL"), "MAXIMUM", "DATA")
    focal_max_large = FocalStatistics(chm, NbrCircle(treetop_radius_large, "CELL"), "MAXIMUM", "DATA")
    # Create a binary raster where 1 = local maxima (treetop) and NoData = everything else. 
    tops = Con(
        ((chm == focal_max_small) | (chm == focal_max_large)) & (chm >= min_canopy_height), 1
    )
    tops.save(treetop_raster)
    # Reduce each local-maxima blob to a single centroid so each tree gets exactly one unique seed value for the watershed.
    tops_polys = os.path.join(scratch_workspace, "tops_polys")
    tops_centroids = os.path.join(scratch_workspace, "tops_centroids")
    arcpy.conversion.RasterToPolygon(tops, tops_polys, "NO_SIMPLIFY")
    arcpy.management.FeatureToPoint(tops_polys, tops_centroids, "CENTROID")
    arcpy.management.AddField(tops_centroids, "SEED_ID", "LONG")
    # Use OBJECTID as a stable unique seed label for allocation.
    arcpy.management.CalculateField(tops_centroids, "SEED_ID", "!OBJECTID!", "PYTHON3")
    # Convert labeled treetop points to a seed raster for EucAllocation.
    arcpy.conversion.PointToRaster(
        tops_centroids, "SEED_ID", treetop_seeds, "MAXIMUM", "", desc.meanCellWidth
    )

    # 3. Allocate canopy pixels to nearest treetop seed
    # EucAllocation is used instead of watershed because broad flat-topped crowns produce near-zero gradients on the inverted CHM
    # 1. Every cell in the analysis space is assigned to the closest treetop seed
    # 2. The output cell value is the value of that nearest SEED_ID
    arcpy.AddMessage("Allocating canopy to nearest treetop...")
    crown_alloc = EucAllocation(Raster(treetop_seeds))
    # Restrict allocation to canopy only (non-canopy remains NoData).
    crown_raster = Con(Raster(canopy_mask), crown_alloc)

    # 4. Convert crown segments to polygons
    arcpy.AddMessage("Converting crown segments to polygons...")
    arcpy.conversion.RasterToPolygon(
        in_raster=crown_raster,
        out_polygon_features=crown_polygons_multi,
        simplify="SIMPLIFY",
        raster_field="VALUE"
    )

    # 5. Multipart to singlepart
    arcpy.AddMessage("Converting multipart polygons to singlepart...")
    arcpy.management.MultipartToSinglepart(crown_polygons_multi, crown_polygons)

    # 6. Add geometry area field if needed and remove tiny polygons
    arcpy.AddMessage("Removing tiny polygons...")
    add_field_if_missing(crown_polygons, "Area_m2", "DOUBLE")
    arcpy.management.CalculateGeometryAttributes(
        crown_polygons,
        [["Area_m2", "AREA"]],
        area_unit="SQUARE_METERS"
    )

    lyr = "crown_lyr"
    arcpy.management.MakeFeatureLayer(crown_polygons, lyr)
    # Remove very small polygons that are usually segmentation artifacts.
    arcpy.management.SelectLayerByAttribute(
        lyr,
        "NEW_SELECTION",
        f'"Area_m2" < {min_crown_area}'
    )
    arcpy.management.DeleteFeatures(lyr)
    arcpy.management.SelectLayerByAttribute(lyr, "CLEAR_SELECTION")

    # 7. Add crown metric fields
    arcpy.AddMessage("Calculating crown width metrics...")
    add_field_if_missing(crown_polygons, "CrownID", "LONG")
    add_field_if_missing(crown_polygons, "MaxWidth_m", "DOUBLE")
    add_field_if_missing(crown_polygons, "MinWidth_m", "DOUBLE")
    add_field_if_missing(crown_polygons, "MeanWidth_m", "DOUBLE")
    add_field_if_missing(crown_polygons, "Perim_m", "DOUBLE")

    arcpy.management.CalculateField(
        crown_polygons,
        "CrownID",
        "!OBJECTID!",
        "PYTHON3"
    )
    arcpy.management.CalculateGeometryAttributes(
        crown_polygons,
        [["Perim_m", "PERIMETER_LENGTH"]],
        length_unit="METERS"
    )

    # Compute area and width metrics from each crown geometry and write to fields.
    with arcpy.da.UpdateCursor(
        crown_polygons,
        ["SHAPE@", "Area_m2", "MaxWidth_m", "MinWidth_m", "MeanWidth_m"]
    ) as cursor:
        for row in cursor:
            geom = row[0]
            area, max_w, min_w, mean_w = polygon_width_metrics(geom)
            row[1] = area
            row[2] = max_w
            row[3] = min_w
            row[4] = mean_w
            cursor.updateRow(row)

    # Build summary statistics for quick QA in the Geoprocessing messages panel.
    mean_widths = [
        row[0] for row in arcpy.da.SearchCursor(crown_polygons, ["MeanWidth_m"])
        if row[0] is not None
    ]
    areas = [
        row[0] for row in arcpy.da.SearchCursor(crown_polygons, ["Area_m2"])
        if row[0] is not None
    ]
    arcpy.AddMessage("=" * 50)
    if areas:
        avg_crown_area_m2 = sum(areas) / len(areas)
        arcpy.AddMessage(f"Mean Crown Area: {avg_crown_area_m2:.2f} m²")

    if mean_widths:
        avg_crown_width_m = sum(mean_widths) / len(mean_widths)
        arcpy.AddMessage(f"Mean Crown Width: {avg_crown_width_m:.2f} m")
    arcpy.AddMessage("=" * 50)

    arcpy.CheckInExtension("Spatial")


def calculate_canopy_metrics(total_area_m2, crown_spacing_sum_m2, num_crowns):
    
    # Area calculations
    remaining_area_m2 = total_area_m2 - crown_spacing_sum_m2
    avg_gap_area_m2 = remaining_area_m2 / num_crowns
    
    # Linear conversions for the spacing between crowns
    square_gap_width_m = math.sqrt(avg_gap_area_m2)
    
    # arcpy.AddMessage results to console
    arcpy.AddMessage("=" * 50)
    arcpy.AddMessage(" FORESTRY CANOPY ANALYSIS REPORT")
    arcpy.AddMessage("=" * 50)
    arcpy.AddMessage(f"Tree Count:               {num_crowns:,} trees")
    arcpy.AddMessage(f"Avg. Gap Area per Tree:   {avg_gap_area_m2:.4f} m²")
    arcpy.AddMessage(f"Linear Gap (Square model): {square_gap_width_m:.2f} m")
    arcpy.AddMessage("=" * 50)


if __name__ == "__main__":
    main()

    total_area_m2 = sum(
        row[0] for row in arcpy.da.SearchCursor(study_area, ["SHAPE@AREA"])
        if row[0] is not None
    )

    crown_spacing_sum_m2 = sum(
        row[0] for row in arcpy.da.SearchCursor(crown_polygons, ["Area_m2"])
        if row[0] is not None
    )

    num_crowns = int(arcpy.management.GetCount(crown_polygons)[0])

    calculate_canopy_metrics(total_area_m2, crown_spacing_sum_m2, num_crowns)