# 0.1.0

- Initial version

# 0.2.0

- add GRB update detection functions and logic
- refactor loader code and other non generic
  parts ([#17](https://github.com/OnroerendErfgoed/brdr/issues/17), ([#7](https://github.com/OnroerendErfgoed/brdr/issues/7))
- add functionality to set a maximum feature size for input
  features (([#50](https://github.com/OnroerendErfgoed/brdr/issues/50))
- fix bug reulting from overlapping features in thematic
  layer ([#46](https://github.com/OnroerendErfgoed/brdr/issues/46))

# 0.2.1

- fixed last_version_date in aligner.get_formula()
- fixed logic of evaluate() in grb.py
- added function to transform geojson to consistent geometry-type (MultiPolygon)