
# Spreadsheet Driven Version
# Date             : Jul 16 2025
#
# Author:          : Jennifer Marr (jennifer.marr@gov.bc.ca)
#
# Purpose:         : Assign crowns to VRI polygons using Identity, then summarize
#                   dominant-overlap species metrics for crown width, count, and spacing.
#
# Inputs:          : crown feature class (must include MeanWidth_m and Area_m2)
#                  : output geodatabase
#                  : VRI polygon feature class (must include LINE_3_TREE_SPECIES)
#                  : optional Excel output path (.xlsx)
#
# Outputs:         : CrownMetrics_By_Species feature class (dissolved by LINE_3_TREE_SPECIES)
#                  : Excel table with species metrics:
#                    LINE_3_TREE_SPECIES, NumCrowns, MeanCrownWidth_m,
#                    CrownSpacing_m, SpeciesArea_m2
#
# History:         : (2026-06-05) Created


"""
Summarize crown metrics by VRI species.

Workflow:
- Run Identity between crowns and VRI polygons.
- Assign each crown to one species by dominant overlap area.
- Summarize NumCrowns, MeanCrownWidth_m, and CrownSpacing_m by LINE_3_TREE_SPECIES.
- Export results to both a feature class and an Excel file.

Requires:
- ArcGIS Pro

Inputs:
- Crown feature class with MeanWidth_m and Area_m2
- VRI polygon feature class with LINE_3_TREE_SPECIES

Outputs:
- Dissolved species feature class with crown metrics
- Excel table of species metrics
"""

import arcpy
import os
import math

arcpy.management.ClearWorkspaceCache() #clean personal arcpro temp folder
arcpy.env.overwriteOutput = True

# -------------------------------------------------------------------
# USER SETTINGS
# -------------------------------------------------------------------
crown_featureclass = arcpy.GetParameterAsText(0)
out_gdb = arcpy.GetParameterAsText(1)
vri_fc = arcpy.GetParameterAsText(2)  
excel_output = arcpy.GetParameterAsText(3)

scratch_workspace = "memory"


def summarize_crowns_by_vri_species(crown_fc, vri_fc, out_fc):
    """
    Identify crowns to VRI polygons and summarize crown metrics by LINE_3_TREE_SPECIES.

    Outputs a dissolved VRI feature class (one polygon per species) with:
    - NumCrowns
    - MeanCrownWidth_m
    - CrownSpacing_m
    """
    required_field = "SPECIES_CD_1"
    vri_fields = [f.name for f in arcpy.ListFields(vri_fc)]
    crown_fields = [f.name for f in arcpy.ListFields(crown_fc)]

    if required_field not in vri_fields:
        raise arcpy.ExecuteError(f"Field {required_field} not found in VRI feature class.")
    if "MeanWidth_m" not in crown_fields or "Area_m2" not in crown_fields:
        raise arcpy.ExecuteError("Crown feature class must contain MeanWidth_m and Area_m2 fields.")

    temp_join = os.path.join("memory", "crowns_vri_identity")
    crown_oid = arcpy.Describe(crown_fc).OIDFieldName

    arcpy.AddMessage("Identifying crowns to VRI polygons using Identity...")
    arcpy.analysis.Identity(
        in_features=crown_fc,
        identity_features=vri_fc,
        out_feature_class=temp_join,
        join_attributes="ALL"
    )

    # Build output geometry by species from identity result (avoiding full VRI processing).
    arcpy.AddMessage("Building VRI species output feature class...")
    arcpy.management.Dissolve(temp_join, out_fc, [required_field])

    add_field_if_missing(out_fc, "NumCrowns", "LONG")
    add_field_if_missing(out_fc, "MeanCrownWidth_m", "DOUBLE")
    add_field_if_missing(out_fc, "CrownSpacing_m", "DOUBLE")

    # Store per-crown metrics so each crown is counted once.
    crown_metrics = {}
    with arcpy.da.SearchCursor(crown_fc, [crown_oid, "MeanWidth_m", "Area_m2"]) as cursor:
        for crownoid, mean_w, area_m2 in cursor:
            crown_metrics[crownoid] = {
                "mean_w": mean_w,
                "area_m2": area_m2 if area_m2 is not None else 0.0,
            }

    # Dominant overlap assignment: pick species with the largest identity piece area per crown.
    crown_to_species = {}
    with arcpy.da.SearchCursor(temp_join, [crown_oid, required_field, "SHAPE@AREA"]) as cursor:
        for crownoid, species, overlap_area in cursor:
            key = species if species not in (None, "") else "Unknown"
            if overlap_area is None:
                overlap_area = 0.0

            if crownoid not in crown_to_species:
                crown_to_species[crownoid] = {
                    "species": key,
                    "overlap": overlap_area,
                }
                continue

            prev = crown_to_species[crownoid]
            if overlap_area > prev["overlap"]:
                crown_to_species[crownoid] = {"species": key, "overlap": overlap_area}
            elif overlap_area == prev["overlap"] and key < prev["species"]:
                # Deterministic tie-breaker
                crown_to_species[crownoid] = {"species": key, "overlap": overlap_area}

    # Summarize by assigned species.
    species_stats = {}
    for crownoid, assigned in crown_to_species.items():
        species = assigned["species"]
        metric = crown_metrics.get(crownoid, {"mean_w": None, "area_m2": 0.0})

        if species not in species_stats:
            species_stats[species] = {
                "num_crowns": 0,
                "sum_mean_width": 0.0,
                "count_mean_width": 0,
                "sum_crown_area": 0.0,
            }

        species_stats[species]["num_crowns"] += 1
        if metric["mean_w"] is not None:
            species_stats[species]["sum_mean_width"] += metric["mean_w"]
            species_stats[species]["count_mean_width"] += 1
        species_stats[species]["sum_crown_area"] += metric["area_m2"]

    # Update output polygons with species summaries.
    excel_rows = []
    with arcpy.da.UpdateCursor(out_fc, [required_field, "SHAPE@AREA", "NumCrowns", "MeanCrownWidth_m", "CrownSpacing_m"]) as cursor:
        for row in cursor:
            species = row[0] if row[0] not in (None, "") else "Unknown"
            species_area_m2 = row[1] if row[1] is not None else 0.0
            stats = species_stats.get(species, None)

            if not stats or stats["num_crowns"] == 0:
                row[2] = 0
                row[3] = None
                row[4] = None
                cursor.updateRow(row)
                excel_rows.append((species, 0, None, None, species_area_m2))
                continue

            num_crowns = stats["num_crowns"]
            mean_crown_width = None
            if stats["count_mean_width"] > 0:
                mean_crown_width = stats["sum_mean_width"] / stats["count_mean_width"]

            remaining_area_m2 = species_area_m2 - stats["sum_crown_area"]
            if remaining_area_m2 < 0:
                remaining_area_m2 = 0.0

            avg_gap_area_m2 = remaining_area_m2 / num_crowns if num_crowns > 0 else 0.0
            crown_spacing_m = math.sqrt(avg_gap_area_m2)

            row[2] = num_crowns
            row[3] = mean_crown_width
            row[4] = crown_spacing_m
            cursor.updateRow(row)
            excel_rows.append((species, num_crowns, mean_crown_width, crown_spacing_m, species_area_m2))

    if not excel_output:
        excel_path = os.path.join(os.path.dirname(out_gdb), "CrownMetrics_By_Species.xlsx")
    else:
        excel_path = excel_output

    export_species_metrics_to_excel(excel_rows, excel_path)

    arcpy.AddMessage("=" * 50)
    arcpy.AddMessage(f"Created species metrics feature class: {out_fc}")
    arcpy.AddMessage(f"Created Excel output: {excel_path}")
    arcpy.AddMessage("=" * 50)


def export_species_metrics_to_excel(rows, excel_path):
    """Write species metrics rows to a temporary table, then export to Excel."""
    temp_table = os.path.join("memory", "species_metrics_table")
    if arcpy.Exists(temp_table):
        arcpy.management.Delete(temp_table)

    arcpy.management.CreateTable("memory", "species_metrics_table")
    add_field_if_missing(temp_table, "LINE_3_TREE_SPECIES", "TEXT")
    add_field_if_missing(temp_table, "NumCrowns", "LONG")
    add_field_if_missing(temp_table, "MeanCrownWidth_m", "DOUBLE")
    add_field_if_missing(temp_table, "CrownSpacing_m", "DOUBLE")
    add_field_if_missing(temp_table, "SpeciesArea_m2", "DOUBLE")

    with arcpy.da.InsertCursor(
        temp_table,
        ["LINE_3_TREE_SPECIES", "NumCrowns", "MeanCrownWidth_m", "CrownSpacing_m", "SpeciesArea_m2"]
    ) as cursor:
        for row in rows:
            cursor.insertRow(row)

    out_dir = os.path.dirname(excel_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    arcpy.conversion.TableToExcel(temp_table, excel_path)


def add_field_if_missing(fc, field_name, field_type):
    """Add a field to a feature class if it doesn't already exist."""
    existing_fields = [f.name for f in arcpy.ListFields(fc)]
    if field_name not in existing_fields:
        arcpy.AddField_management(fc, field_name, field_type)

summarize_crowns_by_vri_species(crown_featureclass, vri_fc, os.path.join(out_gdb, "CrownMetrics_By_Species"))