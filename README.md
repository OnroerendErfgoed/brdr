# `brdr`

a Python library to assist in realigning (multi-)polygons (OGC Simple Features) to reference borders

<!-- badges: start -->

![PyPI - Version](https://img.shields.io/pypi/v/brdr)

<!-- badges: end -->
## Description

### Intro
`brdr` is a Python package that assists in aligning geometric boundaries to reference boundaries. This is an important task in geographic data management to enhance data quality.
* In the context of geographic data management, it is important to have accurate and consistent boundaries for a variety of applications such as calculating areas, analyzing spatial relationships, and visualizing and querying geographic information.
* When creating geographic data, it is often more efficient to derive boundaries from existing reference data rather than collecting new data in the field.
* `brdr` can be used to align boundaries from new data to reference data, ensuring that the boundaries are accurate and consistent.

### Example
The figure below shows:
* the original thematic geometry (blue), 
* A reference layer (yellow-black). 
* The resulting geometry after alignment with `brdr` (green)

![](docs/figures/example.png) 

### Functionalities
`brdr` provides a variety of side-functionalities to assist in aligning boundaries, including:
* Loading thematic data ((Multi-)Polygons): as a dict, geojson or Web Feature Service (WFS-url)
* Loading reference data ((Multi-)Polygons): as a dict, geojson or Web Feature Service (WFS-url)
* (Flanders-specific) Download reference data from GRB-Flanders 
* Align thematic boundaries to reference boundaries
* Calculating a descriptive formulation of a thematic boundary based on a reference layer

### Possible application fields
* Geodata-management:
  * Implementation of `brdr` in business-processes and tooling
  * Bulk geodata-alignment 
  * Alignment after reprojection of data
  * Cleaning data: In a postprocessing-phase, the algorithm executes sliver-cleanup and validity-cleaning on the resulting geometries
  * ... 
* Data-Analysis: Investigate the pattern in deviation and change between thematic and reference boundaries
* Update-detection: Investigate the descriptive formulation before and after alignment to check for (automatic) alignment of geodata 
* ...

### QGIS-script
An implementation of `brdr` for QGIS can be found at [GitHub-brdrQ](https://github.com/OnroerendErfgoed/brdrQ/).
This QGIS- script provides a User Interface to align thematic data to a reference layer, showing the results in the QGIS Table of Contents. 

## Installation

You can install the latest release of `brdr` from
[GitHub](https://github.com/OnroerendErfgoed/brdr/) or
[PyPi](https://pypi.org/project/brdr/):

``` python
pip install brdr
```
### pip-compile
```sh
PIP_COMPILE_ARGS="-v --strip-extras --no-header --resolver=backtracking --no-emit-options --no-emit-find-links"
pip-compile $PIP_COMPILE_ARGS
pip-compile $PIP_COMPILE_ARGS -o requirements-dev.txt --all-extras
```
## Basic example

``` python
from brdr.aligner import Aligner
from shapely import from_wkt
from brdr.enums import OpenbaarDomeinStrategy
from examples import show_results

#CREATE AN ALIGNER
aligner = Aligner(relevant_distance=1, od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
                  threshold_overlap_percentage=50, crs='EPSG:31370')
#ADD A THEMATIC POLYGON TO THEMATIC DICTIONARY
thematic_dict= {"theme_id_1": from_wkt('POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))')}
#ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
reference_dict = {"ref_id_1": from_wkt('POLYGON ((0 1, 0 10,8 10,10 1,0 1))')}
#LOAD THEMATIC DICTIONARY
aligner.load_thematic_data_dict(thematic_dict)
#LOAD REFERENCE DICTIONARY
aligner.load_reference_data_dict(reference_dict)
#EXECUTE THE ALIGNMENT
r, rd, rd_plus, rd_min, sd, si = aligner.process_dict_thematic(relevant_distance=1)
#SHOW RESULTING GEOMETRY AND CHANGES
show_results(r, rd_plus,rd_min,thematic_dict,reference_dict)
```
The resulting figure shows:
* the reference polygon (yellow-black)
* the original geometry (blue)
* the resulting geometry (green line) 
* the added zone (green squares)
* the removed zone (red squares)
<img src="docs/figures/basic_example.png" width="100%" />

More examples can be found in [Examples](https://github.com/OnroerendErfgoed/brdr/tree/main/examples)

## Workflow
(see also Basic example)

To use `brdr`, follow these steps:

* Create a Aligner-class with specific parameters:
  * relevant_distance (m) (default: 1): Distance-parameter used to decide which parts will be aligned, and which parts remain unchanged.
  * od_strategy (enum) (default: SNAP_SINGLE_SIDE): Strategy to align geodata that is not covered by reference-data
  * treshold_overlap_percentage (%)(0-100) (default 50)
  * crs: The Coordinate Reference System (CRS) (default: EPSG:31370 - Belgian Lambert72)
* Load thematic data
* Load reference data
* Process (align) the thematic data
* Results are returned:
  * Resulting geometry
  * Differences: parts that are 'different' from the original geometry (positive or negative)
  * Positive differences: parts that are added to the original geometry
  * Negative differences: parts that are removed form the original geometry
  * Relevant intersections: relevant intersecting parts of the reference geometries
  * Relevant differences: relevant differences of the reference geometries

## The `brdr`-algorithm 

The algorithm for alignment is based on 2 main principles:

* Principle of intentionality: Thematic boundaries can consciously or unconsciously deviate from the reference borders. The algorithm should keep notice of that.
* Selective spatial conservation of shape: The resulting geometry should re-use the shape of the reference borders where aligned is of relevance.

The algorithm can be split into 3 main phases:
* Initialisation: 
  * Deciding which reference polygons are candidate-polygons to re-use its shape. 
  The reference candidate polygons are selected based on spatial intersection with the thematic geometry.
* Processing: 
  * Process all candidate-reference polygons one-by-one
  * Calculate relevant zones for each candidate-reference-polygon
    * relevant intersections: zones that must be present in the final result
    * relevant differences: zones that must be excluded from the final result
![](docs/figures/relevant_zones.png) 
  * Evaluate each candidate based on their relative zones: which parts must be kept and which parts must be excluded
![](docs/figures/evaluate_candidates.png)
  * Union all kept parts to recompose a resulting geometry
* Post-processing: 
  * Technical validation on inner holes and multipolygons
  * Clean-up slivers
  * Make the resulting geometry valid

The figure below shows a schematic overview of the algorithm (in dutch):
![](docs/figures/algorithm.png)

A more in-depth description of the algorithm can be found  in the following article (in dutch):

ADD LINK to report

## Motivation & citation

## Comments and contributions

- Please report any issues or bugs here:
  <https://github.com/OnroerendErfgoed/brdr/issues>.
