# Download-settings: when extracting features by URL
# Limit used when extracting features by URL, using the feature API (f.e. from GRB)
DOWNLOAD_LIMIT = 10000

# value that is used to calculate the boundary of a thematic geometry wherefor the calculation has to be done. (inner part is added)
MAX_OUTER_BUFFER = 50

MAX_SEGMENT_SNAPPING_SIZE = 2  # when real snapping of vertices is used, the input geometry will be split up by default in parts of max X meter
# default CRS:
DEFAULT_CRS = "EPSG:31370"  # BelgianLambert72

# MULTI_SINGLE_ID_SEPARATOR #separator to split multipolygon_ids to single polygons
MULTI_SINGLE_ID_SEPARATOR = "*$*"

PREFIX_FIELDNAME = "brdr_"
BASE_FORMULA_FIELD_NAME = (
    PREFIX_FIELDNAME + "base_formula"
)  # for use in grb_actualisation
FORMULA_FIELD_NAME = PREFIX_FIELDNAME + "formula"
EVALUATION_FIELD_NAME = PREFIX_FIELDNAME + "evaluation"
PREDICTION_SCORE = PREFIX_FIELDNAME + "prediction_score"
PREDICTION_COUNT = PREFIX_FIELDNAME + "prediction_count"
DIFF_PERCENTAGE_FIELD_NAME = PREFIX_FIELDNAME + "diff_percentage"
DIFF_AREA_FIELD_NAME = PREFIX_FIELDNAME + "diff_area"
FULL_BASE_FIELD_NAME = PREFIX_FIELDNAME + "full_base"
FULL_ACTUAL_FIELD_NAME = PREFIX_FIELDNAME + "full_actual"
EQUAL_REFERENCE_FEATURES_FIELD_NAME = PREFIX_FIELDNAME + "equal_reference_features"
OD_ALIKE_FIELD_NAME = PREFIX_FIELDNAME + "od_alike"

AREA_ATTRIBUTE = PREFIX_FIELDNAME + "area"
PERIMETER_ATTRIBUTE = PREFIX_FIELDNAME + "perimeter"
SHAPE_INDEX_ATTRIBUTE = PREFIX_FIELDNAME + "shape_index"

NR_CALCULATION_FIELD_NAME = PREFIX_FIELDNAME + "nr_calculations"
RELEVANT_DISTANCE_FIELD_NAME = PREFIX_FIELDNAME + "relevant_distance"
REMARK_FIELD_NAME = PREFIX_FIELDNAME + "remark"
LAST_VERSION_DATE = "last_version_date"
VERSION_DATE = "version_date"
DATE_FORMAT = "%Y-%m-%d"

# GRB_CONSTANTS
# max buffer (m) around thematic geometry to download reference parcels
GRB_MAX_REFERENCE_BUFFER = 10
# URL of the OGC feature API of actual GRB to extract collections
GRB_FEATURE_URL = "https://geo.api.vlaanderen.be/GRB/ogc/features/collections"
# URL of the OGC feature API of GRB fiscal parcels (situation of 1st of January) to
# extract collections
GRB_FISCAL_PARCELS_URL = "https://geo.api.vlaanderen.be/Adpf/ogc/features/collections"
# Property-name of version_date
GRB_VERSION_DATE = "VERSDATUM"
# Property-name of id of GRB-parcels
GRB_PARCEL_ID = "CAPAKEY"
# Property-name of id of GRB-parcels
GRB_BUILDING_ID = "OIDN"
# Property-name of id of GRB-parcels
GRB_KNW_ID = "OIDN"
