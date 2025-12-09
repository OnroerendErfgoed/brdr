# GRB_CONSTANTS
# max buffer (m) around thematic geometry to download reference parcels
GRB_MAX_REFERENCE_BUFFER = 10
GRB_SUPPORTED_CRS = ["EPSG:31370", "EPSG:3812"]

# URL of the OGC feature API of actual GRB to extract collections
GRB_FEATURE_URL = "https://geo.api.vlaanderen.be/GRB/ogc/features/collections"
# URL of the OGC feature API of GRB fiscal parcels (situation of 1st of January) to
# extract collections
GRB_FISCAL_PARCELS_URL = "https://geo.api.vlaanderen.be/Adpf/ogc/features/collections"
# Property-name of version_date
GRB_VERSION_DATE = "VERSDATUM"
# Property-name of id of GRB-parcels
GRB_PARCEL_ID = "CAPAKEY"
# Property-name of id of GRB-objects
GRB_GENERIC_ID = "OIDN"
