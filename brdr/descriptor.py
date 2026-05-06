from abc import ABC
from abc import abstractmethod
from datetime import datetime
import json
import uuid
from typing import TYPE_CHECKING
from typing import Any

from brdr.constants import MAX_REFERENCE_BUFFER
from brdr.constants import OBSERVATION_FIELD_NAME
from brdr.constants import VERSION_DATE
from brdr.geometry_utils import buffer_pos
from brdr.metadata import get_metadata_observations_from_process_result
from brdr.metadata import reverse_metadata_observations_to_brdr_observation
from brdr.typings import Observation
from brdr.typings import ProcessResult
from brdr.utils import is_brdr_observation
from brdr.utils import urn_from_geom

if TYPE_CHECKING:
    from brdr.aligner import Aligner


class BaseDescriptor(ABC):
    """
    Abstract base class for result description strategies.

    A descriptor enriches a ProcessResult with observations and/or metadata.
    """

    @abstractmethod
    def describe(
        self,
        *,
        aligner: "Aligner",
        thematic_id: Any,
        geometry,
        relevant_distance: float,
        process_result: ProcessResult,
    ) -> ProcessResult:
        """Enrich process_result with descriptor output."""
        pass

    @abstractmethod
    def get_base_observation(
        self, *, feature_properties: dict, metadata_field: str, cache_key: Any = None
    ) -> Observation | None:
        """Decode stored metadata/properties to a BRDR observation."""
        pass

    @abstractmethod
    def get_actual_observation(
        self, *, aligner: "Aligner", process_result: ProcessResult, cache_key: Any = None
    ) -> Observation | None:
        """Resolve actual BRDR observation for a process result."""
        pass


class ObservationSnapshot:
    """
    Lightweight immutable wrapper around a BRDR observation payload.
    """

    def __init__(self, observation: dict):
        self.value = observation


class BaseDescriptorCodec(ABC):
    """Codec contract for descriptor metadata encoding/decoding."""

    @abstractmethod
    def to_metadata_observations(
        self, *, process_result: ProcessResult, reference_lookup: dict
    ) -> list[dict]:
        pass

    @abstractmethod
    def from_metadata(self, metadata: dict) -> dict:
        pass


class BrdrSosaObservationCodec(BaseDescriptorCodec):
    """
    Default BRDR<->SOSA codec used by AlignerDescriptor.
    """

    def to_metadata_observations(
        self, *, process_result: ProcessResult, reference_lookup: dict
    ) -> list[dict]:
        return get_metadata_observations_from_process_result(
            processResult=process_result,
            reference_lookup=reference_lookup,
        )

    def from_metadata(self, metadata: dict) -> dict:
        return reverse_metadata_observations_to_brdr_observation(metadata)


class AlignerDescriptor(BaseDescriptor):
    """
    Default descriptor implementation that adds BRDR observations and metadata.
    """

    def __init__(self, codec: BaseDescriptorCodec | None = None):
        self.codec = codec if codec is not None else BrdrSosaObservationCodec()
        self._base_observation_cache: dict[Any, Observation | None] = {}
        self._actual_observation_cache: dict[Any, Observation | None] = {}

    def describe(
        self,
        *,
        aligner: "Aligner",
        thematic_id: Any,
        geometry,
        relevant_distance: float,
        process_result: ProcessResult,
    ) -> ProcessResult:
        if aligner.add_observations:
            process_result["observations"] = aligner.compare_to_reference(
                process_result.get("result")
            )
            observation_props = aligner.get_observation_properties(process_result)
            props = process_result["properties"]
            props[OBSERVATION_FIELD_NAME] = process_result["observations"]
            props.update(observation_props)
            process_result["properties"] = props

        if aligner.log_metadata:
            actuation_id = uuid.uuid4()
            processor_id = aligner.processor.processor_id.value
            processor_name = type(aligner.processor).__name__
            reference_data = aligner.reference_data
            reference_intersections_ids = reference_data.items.take(
                reference_data.tree.query(buffer_pos(geometry, MAX_REFERENCE_BUFFER))
            ).tolist()
            reference_geometries = []
            for ref_id in reference_intersections_ids:
                feature = reference_data.features[ref_id]
                feat_dict = {
                    "id": feature.brdr_id,
                    "type": f"geo:{feature.geometry.geom_type}",
                    "version_date": reference_data.source.get(VERSION_DATE, ""),
                    "derived_from": {
                        "id": feature.data_id,
                        "type": "geo:Feature",
                        "source": reference_data.source.get("source_url", ""),
                    },
                }
                reference_geometries.append(feat_dict)

            thematic_feature = aligner.thematic_data.features[thematic_id]
            feature_of_interest_id = thematic_feature.brdr_id
            result_urn = urn_from_geom(process_result["result"])
            actuation_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S%z")
            process_result["metadata"] = {}
            process_result["metadata"]["actuation"] = {
                "id": actuation_id.urn,
                "type": "sosa:Actuation",
                "reference_geometries": reference_geometries,
                "changes": "geo:hasGeometry",
                "sosa:hasFeatureOfInterest": {"id": feature_of_interest_id},
                "result": result_urn,
                "result_time": actuation_time,
                "procedure": {
                    "id": processor_id,
                    "implementedBy": processor_name,
                    "type": "sosa:Procedure",
                    "ssn:hasInput": [
                        {
                            "id": "brdr:relevant_distance",
                            "type": "ssn:Input",
                            "input_value": {
                                "type": "xsd:integer",
                                "value": relevant_distance,
                            },
                        },
                    ],
                },
            }
            if process_result["observations"]:
                ref_lookup = reference_data.reference_lookup
                process_result["metadata"]["observations"] = (
                    self.codec.to_metadata_observations(
                        process_result=process_result, reference_lookup=ref_lookup
                    )
                )

        return process_result

    def get_base_observation(
        self, *, feature_properties: dict, metadata_field: str, cache_key: Any = None
    ) -> Observation | None:
        if cache_key is not None and cache_key in self._base_observation_cache:
            return self._base_observation_cache[cache_key]

        base_observation = None
        try:
            if not isinstance(feature_properties, dict):
                raise ValueError("feature_properties must be a dict")
            payload = feature_properties.get(metadata_field)
            if isinstance(payload, str):
                payload = json.loads(payload)
            if is_brdr_observation(payload):
                base_observation = payload
            elif isinstance(payload, dict):
                decoded = self.codec.from_metadata(payload)
                base_observation = (
                    ObservationSnapshot(decoded).value
                    if is_brdr_observation(decoded)
                    else None
                )
        except (
            ValueError,
            TypeError,
            AttributeError,
            json.JSONDecodeError,
        ):
            base_observation = None

        if cache_key is not None:
            self._base_observation_cache[cache_key] = base_observation
        return base_observation

    def get_actual_observation(
        self, *, aligner: "Aligner", process_result: ProcessResult, cache_key: Any = None
    ) -> Observation | None:
        if process_result is None:
            return None

        observation = process_result.get("observations")
        if is_brdr_observation(observation):
            return observation

        geom_process_result = process_result.get("result")
        if geom_process_result is None or geom_process_result.is_empty:
            process_result["observations"] = None
            return None

        if cache_key is None:
            try:
                cache_key = geom_process_result.wkb
            except (AttributeError, TypeError, ValueError):
                cache_key = None

        if cache_key is not None and cache_key in self._actual_observation_cache:
            process_result["observations"] = self._actual_observation_cache[cache_key]
            return process_result["observations"]

        observation = aligner.compare_to_reference(geom_process_result)
        if not is_brdr_observation(observation):
            process_result["observations"] = None
            if cache_key is not None:
                self._actual_observation_cache[cache_key] = None
            return None
        process_result["observations"] = observation
        if cache_key is not None:
            self._actual_observation_cache[cache_key] = observation
        return observation
