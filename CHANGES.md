# 0.1.0

- Initial version

# 0.2.0

- add GRB update detection functions and logic
- refactor loader code and other non generic
  parts ([#17](https://github.com/OnroerendErfgoed/brdr/issues/17), ([#7](https://github.com/OnroerendErfgoed/brdr/issues/7))
- add functionality to set a maximum feature size for input
  features (([#50](https://github.com/OnroerendErfgoed/brdr/issues/50))
- fix bug resulting from overlapping features in thematic
  layer ([#46](https://github.com/OnroerendErfgoed/brdr/issues/46))

# 0.2.1

- fixed last_version_date in aligner.get_formula()
- fixed logic of evaluate() in grb.py
- added function to transform geojson to consistent geometry-type (MultiPolygon)

# 0.3.0

! Not Backwards compatible !

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
- Bug fixing:
    - adding a safe_equals-function to catch GEOsException bug [#71]

# 0.4.0

! Not Backwards compatible !

- Refactoring:
    - Possibility for parallel processing [#97]
    - Changed Aligner constants to init-settings [#83]
    - Refactored ID-handling so strings,integers,... can be used as ID[#110]
    - Cleaned examples [#100]

- Functionalities:
    - Added evaluation-attributes to evaluate()-function [#99]
    - processing-remarks available in geojson-output [#103]
    - Added warning when input/output changed from polygon/multipolygon [#107]

- Bug fixing:
    - Bugfix on version_date [#96]
    - Bugfix on disappearing features [#105]

# 0.5.0

- Refactoring:
    - Updated dependencies [#127]
    - Improvement of the brdr-algorithm implementation [#119]

- Functionalities:
    - *no new functionalities*

- Bug fixing:
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

- Bug fixing:
    - Bugfix for empty geometry[#140]
    - Bugfix for OD-snapping strategies [#146]

# 0.7.0

- Use tracing (=following lines) when snapping to reference borders [#154]
- adding 'partial_snapping' as parameter, to execute post-snapping on geometries that are partially covered [#157]
- Improvement of the od_strategies (SNAP_PREFER_VERTICES, SNAP_NO_PREFERENCE, SNAP_ONLY_VERTICES): performance
  gain [#153]
- Improvement of the od_strategy 'EXCLUDE':performance gain [#150]

# 0.8.0

! Not Backwards compatible !

- Always '0'(zero) added to the relevant distances when calculating predictions, so a (cleaned-up) original geometry is
  available in the process results[#167]
- evaluate()-function: improved and adapted parameters[#164]
- prediction_score: adapted that is an absolute value, expressing the stability-length (in cm)[#162]
- augmented prediction_score when this is stable at max_relevant_distance [#163]
- update_to_actual_grb: added logic that extra parameters (grb_type, max_predictions and full_strategy) can be
  defined[#161]
- Bugfix for updating to actual GRB (max_relevant_distance was always set to the default) [#160]

# 0.8.1

- upgrade hatchling
- added parameter to update_to_actual_grb()

# 0.9.0

- Adapted prediction score to a relative value between 0 and 100 (or -1 = no prediction)[#170]
- bugfix - catch exception when the content of a formula_field cannot be loaded/interpreted[#178]
- Change the default name 'theme_identifier' to 'brdr_id'[#179]
- Added a docker-service to repo - experimental [#143]
- Possibility to align (multi-) linestrings and (multi-)points to reference borders. experimental[#49]

# 0.10.0

! Not Backwards compatible !

- Adapted Enum OpenbaarDomeinStrategy to OpenDomainStrategy
- Adapted Enum Full to FullStrategy
- Added possibilities to align all types of input- and reference-geometries ((multi-)point,line and polygons) [#191]
- Added parameter 'preserve_topology' to keep relations between input-geometries [#90]
- Bugfix on _evaluate-function [#185]

# 0.11.0

! Not Backwards compatible !

- adapted function-name geojson_to_multi
- adapted predictor so point_snapping returns also predictions
- Bugfix for SNAP_ALL_SIDE with big relevant distance (100m) resulting in wrong result [#199]
- prediction-strategy: BEST/ALL/ORIGINAL: check if all records are found in result when using in combination with '
  ONLY_FULL-strategy' [#200]
- GRBActualLoader: support for all types of GRB collections (also lines and points)
- Fix for a better prediction of snapped lines

# 0.12.0

! Not Backwards compatible !

- fix for snapping geometries to points
- fix for snapping geometries to lines
- Changed OpenDomainStrategy from IntEnum to Enum [#210]
- Added enum PredictionStrategy
- New implementation to align points and lines (faster)

# 0.13.0

- removed python 3.9 support [#217]
- Added OSMLoader [#212]
- Added logic to align by lines (process by graph/network) - experimental
- Bugfix: Problem calculating graph_by_multilinestring [#228]