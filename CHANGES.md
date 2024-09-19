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

# 0.3.0

- Refactoring:
  - refactor the structure of the (internal) dicts: dict_series, dict_predicted,.... enhancement question [#57]
  - refactoring of 'formula-unction': generic vs GRB enhancement question [#59]
  - remove deprecated loaders from code enhancement [#77]
  - simplify the core-functionalities of Aligner [#89]
  - cleanup unused functions cleanup [#81]
  - 
- Functionalities:
  - Add brdr-version to formula enhancement [#66]
  - predictor: when multiple predictions with the same resulting geometries are found, only keep the smallest enhancement [#70]
  - Predictor: add a field/attribute that states how many predictions are found, so it can used in the output enhancement [#69]
  - Add function to grb.py: "update_to_actual_version" brdrQ enhancement [#64]
  - GRBSpecificdateLoader: Alignment on GRB (parcels) on specific date [f.e alignment-date] enhancement [#63]
  - Create OnroerendErfgoedLoader; to replace other utirl function to load OE-data enhancement [#65]
- Bugfixing:
  - adding a safe_equals-function to catch GEOsException bug enhancement [#71]
  - research: evaluation of case - to check bug enhancement research [#67]








