from typing import List, Optional, Dict, Any

import numpy as np
import uvicorn
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, model_validator
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType, OpenbaarDomeinStrategy
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader

port = 7999
host = "0.0.0.0"

app = FastAPI()





class ReferenceSource(BaseModel):
    source: str
    version_date: str

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "source": "Adpf",
                    "version_date": "2022-01-01"
                }
            ]
        }
    }

class ReferenceFeature(BaseModel):
    full: bool
    area: float
    percentage: Optional[float]
    version_date: str

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "full": True,
                    "area": 1344.81,
                    "percentage": 100,
                    "version_date": "2019-07-25"
                }
            ]
        }
    }

class Metadata(BaseModel):
    alignment_date: str
    brdr_version: str
    reference_source: ReferenceSource
    full: bool
    area: float
    reference_features: Optional[Dict[Any, ReferenceFeature]]
    #reference_features: Optional[Dict[str, ReferenceFeature]]
    reference_od: Optional[Dict[Any, float]]
    last_version_date: Optional[str]

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "alignment_date": "2025-02-13",
                    "brdr_version": "0.8.1",
                    "reference_source": {
                        "source": "Adpf",
                        "version_date": "2022-01-01"
                    },
                    "full": True,
                    "area": 1344.81,
                    "reference_features": {
                        "24126B0031/00N005": {
                            "full": True,
                            "area": 1344.81,
                            "percentage": 100,
                            "version_date": "2019-07-25"
                        }
                    },
                    "reference_od": None,
                    "last_version_date": "2019-07-25"
                }
            ]
        }
    }

class Properties(BaseModel):
    id: Any
    metadata: Optional[Metadata] = None

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "id": "300",
                    "metadata": {
                        "alignment_date": "2025-02-13",
                        "brdr_version": "0.8.1",
                        "reference_source": {
                            "source": "Adpf",
                            "version_date": "2022-01-01"
                        },
                        "full": True,
                        "area": 1344.81,
                        "reference_features": {
                            "24126B0031/00N005": {
                                "full": True,
                                "area": 1344.81,
                                "percentage": 100,
                                "version_date": "2019-07-25"
                            }
                        },
                        "reference_od": None,
                        "last_version_date": "2019-07-25"
                    }
                }
            ]
        }
    }

class Geometry(BaseModel):
    type: str
    coordinates: List[Any]

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "type": "MultiPolygon",
                    "coordinates": [
                        [
                            [
                                [174111.5042, 179153.9243],
                                [174110.0614, 179154.1094],
                                [174068.867, 179159.3947],
                                [174068.8661, 179159.4262],
                                [174068.8626, 179159.5573],
                                [174073.7483, 179188.9357],
                                [174120.4387, 179180.3235],
                                [174116.1333, 179157.2025],
                                [174111.549, 179153.956],
                                [174111.5042, 179153.9243]
                            ]
                        ]
                    ]
                }
            ]
        }
    }

class Feature(BaseModel):
    type: str
    properties: Properties
    geometry: Geometry

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "type": "Feature",
                    "properties": {
                        "id": "300",
                        "metadata": {
                            "alignment_date": "2025-02-13",
                            "brdr_version": "0.8.1",
                            "reference_source": {
                                "source": "Adpf",
                                "version_date": "2022-01-01"
                            },
                            "full": True,
                            "area": 1344.81,
                            "reference_features": {
                                "24126B0031/00N005": {
                                    "full": True,
                                    "area": 1344.81,
                                    "percentage": 100,
                                    "version_date": "2019-07-25"
                                }
                            },
                            "reference_od": None,
                            "last_version_date": "2019-07-25"
                        }
                    },
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            [
                                [
                                    [174111.5042, 179153.9243],
                                    [174110.0614, 179154.1094],
                                    [174068.867, 179159.3947],
                                    [174068.8661, 179159.4262],
                                    [174068.8626, 179159.5573],
                                    [174073.7483, 179188.9357],
                                    [174120.4387, 179180.3235],
                                    [174116.1333, 179157.2025],
                                    [174111.549, 179153.956],
                                    [174111.5042, 179153.9243]
                                ]
                            ]
                        ]
                    }
                }
            ]
        }
    }

class Params(BaseModel):
    crs: Optional[str] = None
    grb_type: Optional[GRBType] = None
    prediction_strategy: Optional[str] = None

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "crs": "EPSG:31370",
                    "grb_type": "ADP",
                    "prediction_strategy": "ADP"
                }
            ]
        }
    }

class RequestBody(BaseModel):
    features: List[Feature]
    params: Optional[Params] = None

    @model_validator(mode='before')
    def check_unique_ids(cls, values):
        features = values.get('features', [])
        ids = [feature['properties']['id'] for feature in features]
        if len(ids) != len(set(ids)):
            raise ValueError('All feature IDs must be unique')
        return values

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "features": [
                        {
                            "type": "Feature",
                            "properties": {
                                "id": "3"
                            },
                            "geometry":
                                {"coordinates": [
                                    [[174179.0363610595, 179442.4164414695], [174201.26076084137, 179435.4316301095],
                                     [174179.98883533585, 179373.5208021457], [174158.08192697944, 179379.87063065483],
                                     [174179.0363610595, 179442.4164414695]]], "type": "Polygon"}
                        }
                    ],
                    "params": {}
                }

            ]
        }
    }


@app.post("/actualiser")
def actualiser(request_body: RequestBody):
    """
    Create an item with all the information:

    - **name**: each item must have a name
    - **description**: a long description
    - **price**: required
    - **tax**: optional
    """
    try:

        #DEFAULT
        relevant_distances = [
            round(k, 2)
            for k in np.arange(0, 310, 10, dtype=int) / 100
        ]
        crs="EPSG:31370"
        threshold_overlap_percentage=50
        od_strategy=OpenbaarDomeinStrategy.SNAP_ALL_SIDE
        area_limit=100000
        grb_type = GRBType.ADP
        #get geometry
        data_dict ={}
        for f in request_body.features:
            data_dict[f.properties.id] =shape(f.geometry.model_dump())


        #start a new aligner
        aligner = Aligner(crs=crs,
                          threshold_overlap_percentage=threshold_overlap_percentage,
                          od_strategy=od_strategy,
                          area_limit=area_limit
                          )
        #load geometry into thematic dictionary
        aligner.load_thematic_data(DictLoader(
        data_dict = data_dict
        ))
        #load reference data
        aligner.load_reference_data(
            GRBActualLoader(grb_type=grb_type, partition=1000, aligner=aligner)
        )

        #EXECUTE EVALUATION
        aligner.evaluate(relevant_distances=relevant_distances)
        fc = aligner.get_results_as_geojson(
            resulttype=AlignerResultType.EVALUATED_PREDICTIONS,
            formula=True,
            attributes=False,
        )
        #return the featureclasses with the evaluated predictions
        return fc
    except Exception as e:
        raise HTTPException(status_code=400, detail="Missing required fields")

@app.get("/")
def home():
    return ("Welcome to GRB-actualiser webservice!.You can actualise on '/actualiser'."
            "Docs can be found at '/docs'")

def start_server():
    uvicorn.run(app, host=host, port=port)
    print ("Webservice started at " + "http://" + host+":"+ str(port))

if __name__ == "__main__":
    start_server()





