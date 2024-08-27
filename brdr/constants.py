# Thresholds
# Area in m² for excluding candidate reference when overlap(m²) is smaller than the
# threshold
THRESHOLD_EXCLUSION_AREA = 0
# Percentage for excluding candidate reference when overlap(%) is smaller than the
# threshold
THRESHOLD_EXCLUSION_PERCENTAGE = 0

# Buffer parameters:
#   Explanation and examples:
#   https://shapely.readthedocs.io/en/stable/reference/shapely.buffer.html
#   https://postgis.net/docs/ST_Buffer.html
# Distance to limit a buffered corner (MITER-join-style parameter)
MITRE_LIMIT = 10
# Used in buffer-operations to define a quarter circle
QUAD_SEGMENTS = 5

# Correction-parameters (technical)

# Multiplication-factor used in OD-strategy 2 (SNAP-BOTH SIDED) when calculating
# OD-area to take into account
BUFFER_MULTIPLICATION_FACTOR = 1.01
# Threshold-value to exclude circles getting processed (perfect circle = 1) based on
# POLSPY-POPPER algorithm
THRESHOLD_CIRCLE_RATIO = 0.98
# Distance used in a pos_neg_buffer to remove slivers (technical correction)
CORR_DISTANCE = 0.01

# Download-settings: when extracting features by URL
# max buffer around thematic geometry to download reference parcels
MAX_REFERENCE_BUFFER = 10
# Limit used when extracting features by URL, using the feature API (f.e. from GRB)
DOWNLOAD_LIMIT = 10000

# default CRS:
DEFAULT_CRS = "EPSG:31370"  # BelgianLambert72

# MULTI_SINGLE_ID_SEPARATOR #separator to split multipolygon_ids to single polygons
MULTI_SINGLE_ID_SEPARATOR = "*$*"

# URL of the OGC feature API of actual GRB to extract collections
GRB_FEATURE_URL = "https://geo.api.vlaanderen.be/GRB/ogc/features/collections"

# URL of the OGC feature API of GRB fiscal parcels (situation of 1st of January) to extract collections
GRB_FISCAL_PARCELS_URL = "https://geo.api.vlaanderen.be/Adpf/ogc/features/collections"

# input geometry area limit in the selected crs unit (default meters), to prevent resource-intensive operations
AREA_LIMIT = 100000
