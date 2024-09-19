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

! Not Backwards compatable !

- Refactoring:
    - refactor the structure of the (internal) dicts: dict_series, dict_predicted. More logical and faster [#57]
    - refactoring of 'formula-function': more generic [#59]
    - removed deprecated loaders from codebase [#77]
    - simplify the core-functionalities of Aligner: process, predict, compare [#89]
    - cleanup unused functions [#81]

- Functionalities:
    - Add brdr-version to formula [#66]
    - predict: filter duplicate predictions [#70]
    - Predictor: add a attribute (nr_calculations) that states how many predictions are found, so it can used in the
      output [#69]
    - Added GRB-function"update_to_actual_version", to be used in brdrQ [#64]
    - Added GRBSpecificDateLoader: Alignment on GRB (parcels) on specific date [#63]
    - Added OnroerendErfgoedLoader, to load OnroerendErfgoed-data [#65]
- Bugfixing:
    - adding a safe_equals-function to catch GEOsException bug [#71]








