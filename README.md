BRDR
========
BRDR - a Python library to assist in realigning (multi-)polygons (OGC Simple Features) to reference borders

# `brdr`: a Python library to assist in realigning (multi-)polygons (OGC Simple Features) to reference borders

<!-- badges: start -->

<!-- badges: end -->



![](docs/figures/example.png)



## Installation

```sh
PIP_COMPILE_ARGS="-v --strip-extras --no-header --resolver=backtracking --no-emit-options --no-emit-find-links"
pip-compile $PIP_COMPILE_ARGS
pip-compile $PIP_COMPILE_ARGS -o requirements-dev.txt --all-extras
```

You can install the latest release of `brdr` from
[GitHub](https://github.com/OnroerendErfgoed/brdr/) 

``` python
pip install brdr
```


## Basic example

Example to add

``` python
from brdr.aligner import Aligner
from shapely.geometry import Polygon
from examples import show_results

aligner = Aligner()
thematic_dict= {"theme_id_1": Polygon([(0, 0), (0, 9), (5, 10), (10, 0)])}
reference_dict = {"ref_id_1": Polygon([(0, 1), (0, 10), (8, 10), (10, 1)])}
aligner.load_thematic_data_dict(thematic_dict)
aligner.load_reference_data_dict(reference_dict)
r, rd, rd_plus, rd_min, sd, si = aligner.process_dict_thematic(relevant_distance=1)
show_results(r, rd)
```

<img src="docs/figures/example.png" width="100%" />

Explenation: TO ADD

## Getting started

To add

## The algorithm

To add

![](docs/figures/algorithm.png)

## Motivation & citation

## Comments and contributions

- Please report any issues or bugs here:
  <https://github.com/OnroerendErfgoed/brdr/issues>.