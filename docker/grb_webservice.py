import numpy as np
import uvicorn
from fastapi import FastAPI, HTTPException
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType, OpenbaarDomeinStrategy, Full
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
from grb_webservice_typings import ResponseBody, RequestBody

port = 7999
host = "0.0.0.0"

app = FastAPI()


@app.post("/actualiser", response_model=ResponseBody)
def actualiser(request_body: RequestBody):
    """
    Returns GRB-actualised predictions (+ score) for a set of features

    - **featurecollection**: a geojson featurecollection with the features to align
    - **params**: optional: CRS, prediction_strategy
    """
    try:
        # DEFAULT
        relevant_distances = [
            round(k, 2) for k in np.arange(0, 310, 10, dtype=int) / 100
        ]
        crs = "EPSG:31370"
        threshold_overlap_percentage = 50
        od_strategy = OpenbaarDomeinStrategy.SNAP_ALL_SIDE
        area_limit = 100000
        grb_type = GRBType.ADP
        full_strategy = Full.PREFER_FULL
        # get geometry
        data_dict = {}
        for f in request_body.featurecollection.features:
            data_dict[f.id] = shape(f.geometry.model_dump())

        # start a new aligner
        aligner = Aligner(
            crs=crs,
            threshold_overlap_percentage=threshold_overlap_percentage,
            od_strategy=od_strategy,
            area_limit=area_limit,
        )
        # load geometry into thematic dictionary
        aligner.load_thematic_data(DictLoader(data_dict=data_dict))
        # load reference data
        aligner.load_reference_data(
            GRBActualLoader(grb_type=grb_type, partition=1000, aligner=aligner)
        )

        # EXECUTE EVALUATION
        aligner.evaluate(
            relevant_distances=relevant_distances, full_strategy=full_strategy
        )

        # GET RESULTS
        fc = aligner.get_results_as_geojson(
            resulttype=AlignerResultType.EVALUATED_PREDICTIONS,
            formula=True,
            attributes=False,
        )
        # return the featureclasses with the evaluated predictions
        return fc
    except Exception as e:
        raise HTTPException(status_code=400, detail="Missing required fields")


@app.get("/")
def home():
    return (
        "Welcome to GRB-actualiser webservice!.You can actualise on '/actualiser'."
        "Docs can be found at '/docs'"
    )


def start_server():
    uvicorn.run(app, host=host, port=port)
    print("Webservice started at " + "http://" + host + ":" + str(port))


if __name__ == "__main__":
    start_server()
