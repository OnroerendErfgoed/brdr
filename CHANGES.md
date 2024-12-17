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

# 0.4.0

! Not Backwards compatable !

- Refactoring:
    - Possibility for parallel processing [#97]
    - Changed Aligner constants to init-settings [#83]
    - Refactored ID-handling so strings,integers,... can be used as ID[#110]
    - Cleaned examples [#100]

- Functionalities:
    - Added evaluation-attributes to evaluate()-function [#99]
    - processing-remarks available in geojson-output [#103]
    - Added warning when input/output changed from polygon/multipolygon [#107]

- Bugfixing:
    - Bugfix on version_date [#96]
    - Bugfix on disappearing features [#105]

# 0.5.0

- Refactoring:
    - Updated dependencies [#127]
    - Improvement of the brdr-algorithm implementation [#119]

- Functionalities:
    - *no new functionalities*

- Bugfixing:
    - Bugfix needed when using non-textual identifiers (f.e. integers) when searching for GRB changes [#115]
    - Bugfix in predicted geometries: When no stable geometries found an empty list has to be returned instead of
      0.0 [#118]
    - Bugfix when reloading new thematic dictionary; the thematic union has to be recalculated [#122]
    - Bugfix for default list of relevant distances; Range-calculation omits the upper value [#117]
    - Bugfix when using 0m as relevant distance for GRB-actualisation [#121]
    - Bugfix so GRB-actualisation also works with local layers [#120]

# 0.6.0

- Refactoring:
    - Recalculation of the shape_index parameter [#137]
    - Cleaned up several OD-strategies that are not used[#136]
    - Performance gain for the algorithm by only calculating the 'outer boundary' [#142]

- Functionalities:
    - Predictor: Added sorting of predictions based on 'stability length', so the most logical predictions can be
      proposed [#133]
    - Evaluate: aligned features based on 'prediction_score' [#135][#58]
    - Evaluate: Added possibility to prefer 'full results' when evaluating alignments [#144]
    - Added new OD-strategies based on new snapping-function: SNAP_PREFER_VERTICES,
      SNAP_NO_PREFERENCE,SNAP_ONLY_VERTICES [#147]

- Bugfixing:
    - Bugfix for empty geometry[#140]
    - Bugfix for OD-snapping strategies [#146]




