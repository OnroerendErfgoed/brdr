# define a typeddict thematic_data with keys name: str and geom: geometry
from typing import Dict, Any
from typing import List
from typing import TypedDict
from typing import Union

from shapely.geometry.base import BaseGeometry


class GeoJSONGeometry(TypedDict):
    type: str
    coordinates: List


class Crs(TypedDict):
    type: str
    properties: dict


class Feature(TypedDict):
    type: str
    geometry: GeoJSONGeometry
    properties: dict


class FeatureCollection(TypedDict, total=False):
    type: str
    name: str
    crs: Crs
    features: List[Feature]
    __extra_items__: Dict[str, str]


class ProcessResult(TypedDict, total=False):
    """
    A dictionary structure representing the output of a geometric alignment process.

    Attributes
    ----------
    result : BaseGeometry
        The final aligned geometry.
    result_diff : BaseGeometry
        The total geometric difference.
    result_diff_plus : BaseGeometry
        The added geometric area.
    result_diff_min : BaseGeometry
        The removed geometric area.
    properties : Dict[str, Any]
        Calculated metrics and feature attributes.
    metadata : Dict[Any, Any]
        Dictionary containing execution metadata.
    observation : Dict[Any, Any]
        Dictionary containing the alignment observation details.

    Notes
    -----
    The following diagram shows the relationship between the components:

    ```{mermaid}
    graph LR
        A[Original] --> D{Aligner}
        B[Target] --> D
        D --> R[result]
        R --> DP[result_diff_plus]
        R --> DM[result_diff_min]
    ```

    Examples
    --------
    >>> result: ProcessResult = {
    ...     "result": Point(0, 0),
    ...     "properties": {"brdr_stability": True}
    ... }
    """

    result: BaseGeometry
    result_diff: BaseGeometry
    result_diff_plus: BaseGeometry
    result_diff_min: BaseGeometry
    result_relevant_intersection: BaseGeometry
    result_relevant_diff: BaseGeometry
    properties: Dict[str, Any]
    metadata: Dict[Any, Any]
    observation: Dict[Any, Any]


InputId = Union[str, int]


class ObservationReference(TypedDict, total=False):
    area: int
    full: bool
    percentage: int


class Observation(TypedDict, total=False):
    alignment_date: str
    brdr_version: str
    reference_source: str
    full: bool
    area: int
    reference_features: Dict[InputId, ObservationReference]
    reference_od: Dict
