import uuid
from datetime import datetime
from typing import Dict
from typing import List

from brdr.typings import ProcessResult
from brdr.utils import is_brdr_observation

OBS_PROP_AREA_OVERLAP = "brdr:area_overlap"
OBS_PROP_AREA_OVERLAP_PERCENTAGE = "brdr:area_overlap_percentage"
OBS_PROP_AREA_OVERLAP_FULL = "brdr:area_overlap_full"
OBS_PROP_AREA = "brdr:area"
OBS_PROP_AREA_OPEN_DOMAIN = "brdr:area_open_domain"
OBS_PROP_LENGTH_OVERLAP = "brdr:length_overlap"
OBS_PROP_LENGTH_OVERLAP_PERCENTAGE = "brdr:length_overlap_percentage"
OBS_PROP_LENGTH = "brdr:length"
OBS_PROP_LENGTH_OPEN_DOMAIN = "brdr:length_open_domain"
OBS_PROP_COUNT_OVERLAP = "brdr:count_overlap"
OBS_PROP_COUNT_OVERLAP_PERCENTAGE = "brdr:count_overlap_percentage"
OBS_PROP_COUNT = "brdr:count"
OBS_PROP_COUNT_OPEN_DOMAIN = "brdr:count_open_domain"
OBS_PROP_DE9IM = "brdr:de9im_intersection_matrix"

PROC_AREA_OVERLAP = "brdr:observation_procedure_area_overlap"
PROC_AREA_OVERLAP_PERCENTAGE = "brdr:observation_procedure_area_overlap_percentage"
PROC_AREA_OVERLAP_FULL = "brdr:observation_procedure_area_overlap_full"
PROC_AREA = "brdr:observation_procedure_area"
PROC_AREA_OPEN_DOMAIN = "brdr:observation_procedure_area_open_domain"
PROC_DE9IM = "brdr:observation_procedure_de9im_intersection_matrix"


def get_metadata_observations_from_process_result(
    processResult: ProcessResult,
    reference_lookup: Dict[any, any],
) -> List[Dict]:
    observation = processResult["observations"]
    actuation_metadata = processResult["metadata"]["actuation"]
    observation_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S%z")
    sensor_uuid = uuid.uuid4()
    result_id = actuation_metadata["result"]
    observation_metadata = {
        "type": "sosa:Observation",
        "has_feature_of_interest": result_id,
        "made_by_sensor": sensor_uuid.urn,
        "result_time": observation_time,
    }

    observations = []
    measure_type = observation.get("measure_type", "area")
    for ref_id, observations_dict in observation["reference_features"].items():
        reference_id = reference_lookup[ref_id]
        feature_measure_type = observations_dict.get("measure_type", measure_type)
        de9im = observations_dict.get("de9im")
        if de9im is not None:
            observations.append(
                {
                    **observation_metadata,
                    "id": uuid.uuid4().urn,
                    "observed_property": OBS_PROP_DE9IM,
                    "has_feature_of_interest": reference_id,
                    "result": {"value": de9im, "type": "string"},
                    "used_procedure": PROC_DE9IM,
                    "used": result_id,
                }
            )
        area = observations_dict.get("area")
        if area is not None:
            observations.append(
                {
                    **observation_metadata,
                    "id": uuid.uuid4().urn,
                    "observed_property": OBS_PROP_AREA_OVERLAP,
                    "has_feature_of_interest": reference_id,
                    "result": {"value": area, "type": "float"},
                    "used_procedure": PROC_AREA_OVERLAP,
                    "used": result_id,
                }
            )
        percentage = observations_dict.get("percentage")
        if percentage is not None:
            observations.append(
                {
                    **observation_metadata,
                    "id": uuid.uuid4().urn,
                    "has_feature_of_interest": reference_id,
                    "observed_property": OBS_PROP_AREA_OVERLAP_PERCENTAGE,
                    "result": {"value": percentage, "type": "float"},
                    "used_procedure": PROC_AREA_OVERLAP_PERCENTAGE,
                    "used": result_id,
                }
            )
        if feature_measure_type != "area":
            measure_value = observations_dict.get(feature_measure_type)
            if measure_value is not None:
                observations.append(
                    {
                        **observation_metadata,
                        "id": uuid.uuid4().urn,
                        "observed_property": f"brdr:{feature_measure_type}_overlap",
                        "has_feature_of_interest": reference_id,
                        "result": {"value": measure_value, "type": "float"},
                        "used_procedure": f"brdr:observation_procedure_{feature_measure_type}_overlap",
                        "used": result_id,
                    }
                )
            if percentage is not None:
                observations.append(
                    {
                        **observation_metadata,
                        "id": uuid.uuid4().urn,
                        "has_feature_of_interest": reference_id,
                        "observed_property": f"brdr:{feature_measure_type}_overlap_percentage",
                        "result": {"value": percentage, "type": "float"},
                        "used_procedure": f"brdr:observation_procedure_{feature_measure_type}_overlap_percentage",
                        "used": result_id,
                    }
                )
    full = observation.get("full")
    if full is not None:
        observations.append(
            {
                **observation_metadata,
                "id": uuid.uuid4().urn,
                "has_feature_of_interest": result_id,
                "observed_property": OBS_PROP_AREA_OVERLAP_FULL,
                "result": {"value": full, "type": "boolean"},
                "used_procedure": PROC_AREA_OVERLAP_FULL,
            }
        )

    area = observation.get("area")
    if area is not None:
        observations.append(
            {
                **observation_metadata,
                "id": uuid.uuid4().urn,
                "has_feature_of_interest": result_id,
                "observed_property": OBS_PROP_AREA,
                "result": {"value": area, "type": "float"},
                "used_procedure": PROC_AREA,
            }
        )
    if measure_type != "area":
        measure_value = observation.get(measure_type)
        if measure_value is not None:
            observations.append(
                {
                    **observation_metadata,
                    "id": uuid.uuid4().urn,
                    "has_feature_of_interest": result_id,
                    "observed_property": f"brdr:{measure_type}",
                    "result": {"value": measure_value, "type": "float"},
                    "used_procedure": f"brdr:observation_procedure_{measure_type}",
                }
            )
    # Primary key is reference_od; keep area_od fallback for backward compatibility.
    reference_od = observation.get("reference_od")
    if reference_od is None:
        reference_od = observation.get("area_od")
    area_od = None
    if isinstance(reference_od, dict):
        area_od = reference_od.get("area")
    if area_od is not None:
        observations.append(
            {
                **observation_metadata,
                "id": uuid.uuid4().urn,
                "has_feature_of_interest": result_id,
                "observed_property": OBS_PROP_AREA_OPEN_DOMAIN,
                "result": {"value": area_od, "type": "float"},
                "used_procedure": PROC_AREA_OPEN_DOMAIN,
            }
        )
    if measure_type != "area" and isinstance(reference_od, dict):
        metric_od = reference_od.get(measure_type)
        if metric_od is not None:
            observations.append(
                {
                    **observation_metadata,
                    "id": uuid.uuid4().urn,
                    "has_feature_of_interest": result_id,
                    "observed_property": f"brdr:{measure_type}_open_domain",
                    "result": {"value": metric_od, "type": "float"},
                    "used_procedure": f"brdr:observation_procedure_{measure_type}_open_domain",
                }
            )

    return observations


def reverse_metadata_observations_to_brdr_observation(metadata: Dict) -> Dict:
    """
    Reconstruct a BRDR observation dictionary from a list of Linked Data observations.
    """
    if not metadata or not "actuation" in metadata or not "observations" in metadata:
        return {}
    observations = metadata["observations"]
    actuation_metadata = metadata["actuation"]
    reverse_ref_lookup = {}
    for ref_geom in actuation_metadata.get("reference_geometries", []):
        uuid_id = ref_geom.get("id")
        original_id = ref_geom.get("derived_from", {}).get("id")
        if uuid_id and original_id:
            reverse_ref_lookup[uuid_id] = original_id

    main_result_id = actuation_metadata.get("result")
    reconstructed = {
        "reference_features": {},
        "full": None,
        "area": None,
        "length": None,
        "count": None,
        "measure_type": "area",
    }
    for obs in observations:
        prop = obs.get("observed_property")
        val = obs.get("result", {}).get("value")
        foi = obs.get("has_feature_of_interest")

        if foi == main_result_id:
            if prop == OBS_PROP_AREA_OVERLAP_FULL:
                reconstructed["full"] = val
            elif prop == OBS_PROP_AREA:
                reconstructed["area"] = val
                reconstructed["measure_type"] = "area"
            elif prop == OBS_PROP_LENGTH:
                reconstructed["length"] = val
                reconstructed["measure_type"] = "length"
            elif prop == OBS_PROP_COUNT:
                reconstructed["count"] = val
                reconstructed["measure_type"] = "count"
            elif prop == OBS_PROP_AREA_OPEN_DOMAIN:
                reconstructed["reference_od"] = {"area": val}
            elif prop == OBS_PROP_LENGTH_OPEN_DOMAIN:
                reconstructed["reference_od"] = {"length": val}
            elif prop == OBS_PROP_COUNT_OPEN_DOMAIN:
                reconstructed["reference_od"] = {"count": val}

        elif foi in reverse_ref_lookup:
            ref_id = reverse_ref_lookup[foi]
            if ref_id not in reconstructed["reference_features"]:
                reconstructed["reference_features"][ref_id] = {}

            if prop == OBS_PROP_DE9IM:
                reconstructed["reference_features"][ref_id]["de9im"] = val
            elif prop == OBS_PROP_AREA_OVERLAP:
                reconstructed["reference_features"][ref_id]["area"] = val
                reconstructed["reference_features"][ref_id]["measure_type"] = "area"
            elif prop == OBS_PROP_AREA_OVERLAP_PERCENTAGE:
                reconstructed["reference_features"][ref_id]["percentage"] = val
                if val == 100:
                    reconstructed["reference_features"][ref_id]["full"] = True
                else:
                    reconstructed["reference_features"][ref_id]["full"] = False
            elif prop == OBS_PROP_LENGTH_OVERLAP:
                reconstructed["reference_features"][ref_id]["length"] = val
                reconstructed["reference_features"][ref_id]["measure_type"] = "length"
            elif prop == OBS_PROP_LENGTH_OVERLAP_PERCENTAGE:
                reconstructed["reference_features"][ref_id]["percentage"] = val
                reconstructed["reference_features"][ref_id]["full"] = val == 100
                reconstructed["reference_features"][ref_id]["measure_type"] = "length"
                reconstructed["measure_type"] = "length"
            elif prop == OBS_PROP_COUNT_OVERLAP:
                reconstructed["reference_features"][ref_id]["count"] = val
                reconstructed["reference_features"][ref_id]["measure_type"] = "count"
            elif prop == OBS_PROP_COUNT_OVERLAP_PERCENTAGE:
                reconstructed["reference_features"][ref_id]["percentage"] = val
                reconstructed["reference_features"][ref_id]["full"] = val == 100
                reconstructed["reference_features"][ref_id]["measure_type"] = "count"
                reconstructed["measure_type"] = "count"

    reconstructed["alignment_date"] = actuation_metadata.get("result_time", None)
    reconstructed["reference_source"] = None
    reconstructed["brdr_version"] = None
    if "reference_od" not in reconstructed:
        reconstructed["reference_od"] = None
    if reconstructed["reference_od"] is not None:
        reconstructed["area_od"] = reconstructed["reference_od"]

    if not is_brdr_observation(reconstructed):
        raise ValueError("The reconstructed object is not a valid brdr observation")

    return reconstructed
