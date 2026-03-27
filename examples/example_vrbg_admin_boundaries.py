import requests

from brdr.aligner import Aligner
from brdr.be.be import BeAdministrativeBoundaryLoader
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader


def _discover_vrbg_collection() -> str:
    """
    Discover a VRBG collection id.

    Preference:
    1. A collection id containing 'gemeente' (municipality-like naming).
    2. The first collection id returned by the API.
    """
    url = "https://geo.api.vlaanderen.be/VRBG/ogc/features/v1/collections"
    response = requests.get(url, params={"f": "json"}, timeout=30)
    response.raise_for_status()
    data = response.json()
    collections = [c.get("id") for c in data.get("collections", []) if c.get("id")]
    if not collections:
        raise ValueError("No VRBG collections found")

    preferred = next((c for c in collections if "gemeente" in c.lower()), None)
    return preferred or collections[0]


if __name__ == "__main__":
    """
    Example: load Belgian administrative boundaries from VRBG as reference data.
    """
    aligner = Aligner(crs="EPSG:3812")

    thematic_id = "my_contour_id"
    thematic_dict = {
        thematic_id: geom_from_wkt(
            "POLYGON ((172450 172250, 172750 172250, 172750 172550, 172450 172550, 172450 172250))"
        )
    }
    aligner.load_thematic_data(DictLoader(data_dict=thematic_dict))

    collection_id = _discover_vrbg_collection()
    print(f"Using VRBG collection: {collection_id}")

    loader = BeAdministrativeBoundaryLoader(
        aligner=aligner,
        collection=collection_id,
        id_property="OIDN",
        partition=1000,
        limit=10000,
        request_timeout=60,
    )
    aligner.load_reference_data(loader)

    relevant_distance = 5
    aligner_result = aligner.process(relevant_distances=[relevant_distance])
    process_results = aligner_result.get_results(aligner=aligner)
    print("result:", process_results[thematic_id][relevant_distance]["result"].wkt)
