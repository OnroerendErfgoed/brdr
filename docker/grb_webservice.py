import numpy as np
import uvicorn
from fastapi import FastAPI, HTTPException
from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType, OpenbaarDomeinStrategy
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
#from brdr.utils import geojson_geometry_to_shapely


port = 7999
host = "0.0.0.0"

app = FastAPI()

from pydantic import BaseModel, Field, model_validator, ValidationError
from typing import List, Dict, Any
from shapely.geometry import shape

class GeoJSONPolygon(BaseModel):
    type: str = Field(..., example="Polygon")
    coordinates: List[List[List[float]]] = Field(
        ...,
        example=[
            [
                [173933.56907947885, 179488.30342623874],
                [173936.35215775482, 179488.67450334204],
                [173934.272076098, 179506.26064825893],
                [173933.01475345079, 179506.0509856817],
                [173931.93712825546, 179505.87128823],
                [173930.62761179, 179520.14188669],
                [173930.46313179, 179521.93427069],
                [173930.39100379, 179522.72044669],
                [173938.4672918, 179523.72991869],
                [173945.75362781016, 179524.6406386901],
                [173945.87503581008, 179523.31763069],
                [173946.02997980008, 179521.62924669],
                [173947.10466781, 179509.9178866799],
                [173947.56194781009, 179504.93478267992],
                [173951.4300438099, 179462.78214265],
                [173944.1195798101, 179462.13452665],
                [173936.12290780008, 179461.42611064983],
                [173935.80130780008, 179464.81087865008],
                [173933.56907947885, 179488.30342623874]
            ]
        ]
    )

class Params(BaseModel):
    crs: str = Field(..., example="EPSG:31370")
    grb_type: str = Field(..., example="ADP")

class RequestBody(BaseModel):
    geometry: GeoJSONPolygon
    params: Params

    @model_validator(mode='before')
    def check_polygon_area(cls, values: Dict[str, Any]):
        geometry = values.get('geometry')
        polygon = shape(geometry)
        area = polygon.area
        if area >= 100000:
            raise ValidationError(f"Polygon area is too large: {area} m²")
        return values

    class Config:
        schema_extra = {
            "example": {
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [173933.56907947885, 179488.30342623874],
                            [173936.35215775482, 179488.67450334204],
                            [173934.272076098, 179506.26064825893],
                            [173933.01475345079, 179506.0509856817],
                            [173931.93712825546, 179505.87128823],
                            [173930.62761179, 179520.14188669],
                            [173930.46313179, 179521.93427069],
                            [173930.39100379, 179522.72044669],
                            [173938.4672918, 179523.72991869],
                            [173945.75362781016, 179524.6406386901],
                            [173945.87503581008, 179523.31763069],
                            [173946.02997980008, 179521.62924669],
                            [173947.10466781, 179509.9178866799],
                            [173947.56194781009, 179504.93478267992],
                            [173951.4300438099, 179462.78214265],
                            [173944.1195798101, 179462.13452665],
                            [173936.12290780008, 179461.42611064983],
                            [173935.80130780008, 179464.81087865008],
                            [173933.56907947885, 179488.30342623874]
                        ]
                    ]
                },
                "params": {
                    "crs": "EPSG:31370",
                    "grb_type": "ADP"
                }
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
        geojson_geometry = request_body.geometry.model_dump()

        #start a new aligner
        aligner = Aligner(crs=crs,
                          threshold_overlap_percentage=threshold_overlap_percentage,
                          od_strategy=od_strategy,
                          area_limit=area_limit
                          )
        #load geometry into thematic dictionary
        aligner.load_thematic_data(DictLoader(
            #data_dict={'id_1': geojson_geometry_to_shapely(geojson_geometry)}
        data_dict = {'id_1': shape(geojson_geometry)}
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
    return "Welcome to GRB actualiser webservice!.You can actualise on '/actualiser'"

# @app.get("/items/{item_id}")
# def read_item(item_id: int, q: Union[str, None] = None):
#     return {"item_id": item_id, "q": q}


def start_server():
    uvicorn.run(app, host=host, port=port)
    print ("Webservice started at " + "http://" + host+":"+ str(port))

if __name__ == "__main__":
    start_server()





