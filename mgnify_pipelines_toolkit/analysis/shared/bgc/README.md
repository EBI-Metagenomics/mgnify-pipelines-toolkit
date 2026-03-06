**Biosynthetic Gene Cluster (BGC) integration and annotation toolkit**

## Overview

The `bgc` module provides tools for integrating and standardizing BGC predictions from multiple specialized tools (GECCO, antiSMASH, SanntiS) into a unified GFF3 format with antiSMASH-compatible JSON sideloader output.

### Purpose

- **Integration**: Merge overlapping BGC predictions from multiple tools into consensus regions
- **Standardization**: Produce GFF3 and JSON outputs compatible with downstream analysis pipelines
- **Enrichment**: Annotate CDS features with BGC membership and functional predictions
- **Interoperability**: Generate antiSMASH sideloader JSON for visualization and further analysis

## Architecture

The module is organized into focused, single-responsibility components:

```
bgc/
├── cli.py              # Command-line interface and workflow orchestration
├── models.py           # Core data models (BGCRegion, CDSRec, MergedRegion)
├── tool_parsers.py     # Tool-specific GFF parsers (GECCO, antiSMASH, SanntiS)
├── merge.py            # Region merging and CDS annotation logic
├── gff_output.py       # GFF3 output formatting
├── sideload_json.py    # antiSMASH JSON sideloader generation and validation
├── README.md           # This file
└── antismash_sideload_schemas
    ├── antismash_schema.py # Pydantic models for validation (auto-generated)
    ├── build_models.py     # Script to regenerate models from antiSMASH schemas
    └── README.md           # AntiSMASH sideloader schema models readme file

```

### Module Descriptions

#### `models.py`
Core data structures for representing BGC regions and CDS features:
- **`BGCRegion`**: Single BGC prediction from a tool (contig, coordinates, attributes, tool name)
- **`CDSRec`**: CDS feature with genomic coordinates and original GFF line
- **`MergedRegion`**: Consensus region formed by merging overlapping BGC predictions

#### `tool_parsers.py`
Abstract base class (`BGCToolParser`) and concrete implementations for parsing tool-specific GFF formats:
- **`GECCOParser`**: Extracts BGC features with Type attributes (e.g., "NRP,Polyketide")
- **`AntiSMASHParser`**: Parses region and gene-level annotations (product, gene_functions, etc.)
- **`SanntiSParser`**: Extracts CLUSTER features with MiBIG similarity annotations

Each parser returns:
1. List of `BGCRegion` objects
2. Dictionary of gene-level annotations (antiSMASH only)

#### `merge.py`
Consensus region creation and CDS annotation:
- **`merge_overlaps()`**: Merges overlapping BGC predictions into `MergedRegion` objects
- **`support_and_filter_cds()`**: Assigns CDS features to merged regions and adds BGC-specific attributes

Merging strategy:
- Regions from the same contig that overlap by ≥1 bp are merged
- Merged regions track all contributing member predictions
- CDS features are retained if they fall within merged region boundaries

#### `gff_output.py`
GFF3 output generation:
- **`build_region_lines()`**: Converts `MergedRegion` objects to GFF3-formatted region features
- Sorts output by contig and coordinate for standard compliance
- Adds BGC-specific attributes (bgc_id, bgc_tools, member counts, etc.)

#### `sideload_json.py`
antiSMASH-compatible JSON sideloader generation and validation:
- **`build_sideload_json_payload()`**: Constructs JSON payload from merged regions (subregions-only format)
- **`write_sideload_json()`**: Writes JSON file with optional schema validation
- **`validate_sideload_json()`**: Validates using Pydantic models for fast, type-safe validation

**Key features**:
- GFF coordinates (1-based inclusive) → JSON coordinates (0-based start, end-exclusive)
- Tool metadata (name, version, description, timestamp)
- Per-subregion details (ID, tools, sources, member attributes)
- **Opt-in validation**: Use `--validate_json` flag to enable schema validation (disabled by default for performance)
- **Pydantic-based validation**: Fast, type-safe validation based on `antismash_sideload_schemas/antismash_schema.py` (no external schema files required)

#### `antismash_sideload_schemas/antismash_schema.py`
Auto-generated Pydantic models from official antiSMASH sideloader JSON schema:
- **`Model`**: Root model for validating complete sideloader payloads
- **`Tool`**: Tool metadata validation
- **`Record`**: Record-level validation with subregions/protoclusters
- **`Subregion`**, **`Protocluster`**: BGC region type validators
- **`Details`**: Flexible key-value attribute validation

#### `cli.py`
Command-line interface and workflow orchestration:
- Validates input files
- Loads base CDS from reference GFF
- Parses BGC predictions from optional tool GFFs
- Merges overlapping regions
- Writes integrated GFF3 and JSON outputs

## Usage

### Command-Line Interface

```bash
python -m mgnify_pipelines_toolkit.analysis.shared.bgc.cli \
  --base_gff base_annotations.gff \
  --gecco_gff gecco_output.gff \
  --antismash_gff antismash_output.gff \
  --sanntis_gff sanntis_output.gff \
  --output_gff integrated_bgc.gff \
  [--validate_json] \
  [--verbose INFO]
```

**Required arguments**:
- `--base_gff`: Reference GFF3 with CDS annotations (may be compressed with gzip)
- `--output_gff`: Path for integrated output GFF3

**Optional tool inputs** (at least one required):
- `--gecco_gff`: GECCO BGC predictions
- `--antismash_gff`: antiSMASH BGC predictions
- `--sanntis_gff`: SanntiS BGC predictions

**Flags**:
- `--validate_json`: Enable schema validation for JSON sideloader output (default: disabled)
- `--verbose`: Logging level (DEBUG, INFO, WARNING, ERROR)

### Output Files

The tool produces two output files with matching basenames:

1. **Integrated GFF3** (`output.gff`)
   - Base CDS features filtered to BGC regions
   - BGC region features with merged attributes
   - CDS features annotated with BGC membership and tool-specific attributes
   - Sorted by contig and coordinate

2. **antiSMASH Sideloader JSON** (`output.json`)
   - Automatically derived from output GFF path (e.g., `output.gff` → `output.json`)
   - Compatible with antiSMASH visualization tools
   - Contains subregions with detailed annotations
   - Optionally validated against antiSMASH schema using Pydantic

### Programmatic Usage

```python
from pathlib import Path
from mgnify_pipelines_toolkit.analysis.shared.bgc.tool_parsers import GECCOParser, AntiSMASHParser
from mgnify_pipelines_toolkit.analysis.shared.bgc.merge import merge_overlaps, support_and_filter_cds
from mgnify_pipelines_toolkit.analysis.shared.bgc.sideload_json import write_sideload_json

# Parse BGC predictions
gecco_parser = GECCOParser()
antismash_parser = AntiSMASHParser()

gecco_regions, _ = gecco_parser.parse_regions(Path("gecco.gff"))
antismash_regions, gene_ann = antismash_parser.parse_regions(Path("antismash.gff"))

all_regions = gecco_regions + antismash_regions

# Merge overlapping regions
merged_regions = merge_overlaps(all_regions)

# Generate antiSMASH sideloader JSON (with validation)
write_sideload_json(
    out_json=Path("output.json"),
    merged_regions=merged_regions,
    validate=True,  # Enable Pydantic validation
    tool_name="Custom BGC Analyzer",
    tool_version="1.0.0",
    tool_description="Custom BGC integration workflow"
)
```

## Input Requirements

### Base GFF Format
Must contain CDS features with standard GFF3 formatting:
```
##gff-version 3
contig1  Prodigal  CDS  100  500  .  +  0  ID=cds_001;product=hypothetical protein
contig1  Prodigal  CDS  600  900  .  +  0  ID=cds_002;product=oxidoreductase
```

### Tool-Specific GFF Formats

#### GECCO
```
##gff-version 3
contig1  GECCO v0.9.8  BGC  100  1500  .  .  .  ID=cluster1;Type=NRP,Polyketide
```
**Required attributes**: `Type` (BGC classification)

#### antiSMASH
```
##gff-version 3
contig1  antiSMASH:8.0.1  region  100  2000  .  .  .  ID=region1;product=sactipeptide
contig1  antiSMASH:8.0.1  gene    150  450   .  +  .  ID=gene1;Parent=region1;gene_functions=biosynthetic...
```
**Required attributes**: `product` (region), `gene_functions` (gene-level)

#### SanntiS
```
##gff-version 3
contig1  SanntiSv0.9.3.3  CLUSTER  100  1800  .  .  .  nearest_MiBIG=BGC0001234;nearest_MiBIG_class=NRP
```
**Required attributes**: `nearest_MiBIG` (at minimum)

## Merging Strategy

### Region Overlap Rules
Two BGC regions are merged if:
1. They are on the **same contig**
2. Their coordinates **overlap by ≥1 bp**

Example:
```
Region A: [100, 500]  (GECCO)
Region B: [450, 800]  (antiSMASH)
Merged:   [100, 800]  (contains both members)
```

### CDS Filtering
CDS features are **retained** if:
1. They fall **completely within** a merged BGC region's boundaries
2. They are assigned BGC-specific attributes from contributing tools

### Attribute Inheritance
- **Union**: All unique attribute values from member regions are preserved
- **Counts**: Track number of tools and member BGCs
- **Labels**: Prioritize meaningful labels (e.g., GECCO Type, antiSMASH product)

## JSON Sideloader Format

The antiSMASH sideloader JSON follows the official schema:

```json
{
  "tool": {
    "name": "MGnify BGC mapper",
    "version": "mgnify_pipelines_toolkit_v1.4.18",
    "description": "BGC subregions integrated from GECCO/antiSMASH/SanntiS into a base GFF.",
    "configuration": {}
  },
  "records": [
    {
      "name": "contig1",
      "subregions": [
        {
          "start": 99,
          "end": 1500,
          "label": "NRP,Polyketide",
          "details": {
            "ID": "contig1|bgc:100-1500",
            "bgc_tools": ["gecco", "antismash"],
            "member_bgcs": "2",
            "gecco_bgc_type": ["NRP,Polyketide"],
            "antismash_product": ["sactipeptide"]
          }
        }
      ]
    }
  ],
  "timestamp": "2026-03-05T12:48:15.123456+00:00"
}
```

### Schema Validation

Validation is performed using **Pydantic models** auto-generated from the official antiSMASH sideloader schema. This provides fast, type-safe validation with clear error messages.

**Validation behavior**:
- JSON is **always written** first (you get output even if validation fails)
- If `--validate_json` is enabled, Pydantic validation runs after writing
- Validation errors raise `ValueError` with detailed field-level error messages

**Pydantic schema location**: `mgnify_pipelines_toolkit/analysis/shared/bgc/antismash_sideload_schemas/antismash_schema.py`

**Example validation error**:
```
ValueError: Schema validation failed:
2 validation errors for Model
records.0.subregions.0.label
  String should have at most 20 characters [type=string_too_long]
tool.name
  Field required [type=missing]
```

### Legacy JSON Schema Files

The original JSON schema files are still available in `antismash_sideload_schemas/general/` for reference, but are **no longer used** for validation. These can be safely removed after confirming the Pydantic validation works with your data.

## Dependencies

### Required
- Python 3.8+
- `pydantic` >= 2.0: For JSON schema validation (if using `--validate_json`)

### Optional
- `jsonschema` + `referencing`: **No longer required** (replaced by Pydantic validation)

## Testing

Unit tests are located in `tests/unit/analysis/shared/bgc/`:

```bash
# Run all BGC tests
pytest tests/unit/analysis/shared/bgc/ -v

# Run specific test modules
pytest tests/unit/analysis/shared/bgc/test_tool_parsers.py -v
pytest tests/unit/analysis/shared/bgc/test_sideload_json.py -v
```

## Troubleshooting

### Common Issues

#### "At least one optional predictor GFF must be provided"
**Cause**: No tool GFF files were specified
**Solution**: Provide at least one of `--gecco_gff`, `--antismash_gff`, or `--sanntis_gff`

#### "No BGC regions parsed from optional inputs"
**Cause**: Tool GFFs don't contain expected feature types or attributes
**Solution**: Check that:
- GECCO GFFs have `BGC` features with `Type` attributes
- antiSMASH GFFs have `region` features with `product` attributes
- SanntiS GFFs have `CLUSTER` features with `nearest_MiBIG` attributes

#### "No module named 'pydantic'"
**Cause**: `--validate_json` flag used but `pydantic` not installed
**Solution**: Install pydantic:
```bash
pip install pydantic
# OR
# Remove --validate_json flag to skip validation
```

#### JSON validation fails with Pydantic errors
**Cause**: Generated JSON doesn't conform to antiSMASH schema
**Solution**:
1. Check error message for specific field that failed validation
2. Verify that merged regions have valid attributes
3. Ensure tool metadata (name, version) are set correctly
4. Common issues:
   - `label` field exceeds 20 characters (truncated automatically)
   - `records` list is empty (need at least 1 record)
   - Required fields missing from `tool` section

### Debug Mode

Enable detailed logging:
```bash
python -m mgnify_pipelines_toolkit.analysis.shared.bgc.cli \
  --base_gff base.gff \
  --gecco_gff gecco.gff \
  --output_gff output.gff \
  --verbose DEBUG
```

## Design Principles

1. **Separation of Concerns**: Each module has a single, well-defined responsibility
2. **Extensibility**: New tool parsers can be added by implementing `BGCToolParser`
3. **Fail-Safe**: JSON is always written before validation (artifacts even if validation fails)
4. **Performance**: Schema validation is opt-in to avoid unnecessary overhead; Pydantic validation is ~10x faster than JSON schema
5. **Standards Compliance**: Outputs conform to GFF3 and antiSMASH specifications
6. **Type Safety**: Pydantic models provide compile-time type checking and IDE support


## Contributing

When adding support for new BGC prediction tools:

1. **Implement `BGCToolParser`** in `tool_parsers.py`:
   ```python
   class NewToolParser(BGCToolParser):
       @property
       def tool_name(self) -> str:
           return "newtool"

       def parse_regions(self, gff_path: Path) -> tuple[list[BGCRegion], dict[str, dict[str, str]]]:
           # Parse tool-specific GFF format
           # Return (regions, gene_annotations)
   ```

2. **Register parser** in `cli.py`:
   ```python
   _PARSERS: dict[str, BGCToolParser] = {
       "gecco": GECCOParser(),
       "sanntis": SanntiSParser(),
       "antismash": AntiSMASHParser(),
       "newtool": NewToolParser(),  # Add here
   }
   ```

3. **Add CLI option** in `cli.py`:
   ```python
   @click.option("--newtool_gff", type=click.Path(path_type=Path), default=None, help="Optional NewTool output GFF.")
   ```

4. **Write tests** in `tests/unit/analysis/shared/bgc/test_tool_parsers.py`

## License

Copyright 2024-2026 EMBL - European Bioinformatics Institute

Licensed under the Apache License, Version 2.0
