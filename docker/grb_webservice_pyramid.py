from wsgiref.simple_server import make_server

import numpy as np
from pyramid.config import Configurator
from pyramid.response import Response
from pyramid.view import view_config

from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType, OpenbaarDomeinStrategy
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
from brdr.utils import geojson_geometry_to_shapely

host = "0.0.0.0"
host = "localhost"
port = 8765


@view_config(route_name="home")
def home(request):
    return Response(
        "Welcome to GRB actualiser webservice!.You can actualise on '/actualiser'"
    )


@view_config(route_name="actualiser", renderer="json", request_method="POST")
def actualiser(request):
    """

    :param request:
    :return:
    """
    try:
        data = request.json_body
        # DEFAULT
        relevant_distances = [
            round(k, 2) for k in np.arange(0, 310, 10, dtype=int) / 100
        ]
        crs = "EPSG:31370"
        threshold_overlap_percentage = 50
        od_strategy = OpenbaarDomeinStrategy.SNAP_ALL_SIDE
        area_limit = 100000
        grb_type = GRBType.ADP

        # get geometry
        geojson_geometry = data["geometry"]

        # get parameters
        if "params" in data:
            params = data["params"]
            if "crs" in params:
                crs = params["crs"]
            if "relevant_distances" in params:
                relevant_distances = params["relevant_distances"]
            if "od_strategy" in params:
                od_strategy = params["od_strategy"]
            if "grb_type" in params:
                grb_type = GRBType[params["grb_type"]]

        # validate geometry & parameters
        # TODO

        # start a new aligner
        aligner = Aligner(
            crs=crs,
            threshold_overlap_percentage=threshold_overlap_percentage,
            od_strategy=od_strategy,
            area_limit=area_limit,
        )
        # load geometry into thematic dictionary
        aligner.load_thematic_data(
            DictLoader(
                data_dict={"id_1": geojson_geometry_to_shapely(geojson_geometry)}
            )
        )
        # load reference data
        aligner.load_reference_data(
            GRBActualLoader(grb_type=grb_type, partition=1000, aligner=aligner)
        )

        # EXECUTE EVALUATION
        aligner.evaluate(relevant_distances=relevant_distances)
        fc = aligner.get_results_as_geojson(
            resulttype=AlignerResultType.EVALUATED_PREDICTIONS,
            formula=True,
            attributes=False,
        )
        # return the featureclasses with the evaluated predictions
        return fc
    except Exception as e:
        return Response(str(e), status=400)


def start_service():
    with Configurator() as config:
        config.add_route("home", "/")
        config.add_route("actualiser", "/actualiser")
        config.scan()
        app = config.make_wsgi_app()
    server = make_server(host, port, app)
    print("webservice active on http://localhost:" + str(port))
    server.serve_forever()


if __name__ == "__main__":
    start_service()
