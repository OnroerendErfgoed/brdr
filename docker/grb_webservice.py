from wsgiref.simple_server import make_server

import numpy as np
from pyramid.config import Configurator
from pyramid.response import Response
from pyramid.view import view_config
from shapely.geometry.geo import shape

from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader

#from shapely.geometry.geo import shape, mapping

host = "0.0.0.0"
port = 8000

@view_config(
    route_name="home"
)
def home(request):
    return Response("Welcome to GRB actualiser webservice!.You can actualise on '/actualiser'")

@view_config(route_name='actualiser', renderer='json', request_method='POST')
def actualiser(request):
    try:
        #get geometry and parameters from request body
        data = request.json_body
        geometry = data['geometry']
        if 'params' in data:
            params = data['params']
        else:
            params = None

        relevant_distances = params["relevant_distances"] if (params is not None and 'relevant_distances' in params) else [
            round(k, 2)
            for k in np.arange(0, 310, 10, dtype=int) / 100
        ]

        #start a new aligner
        aligner = Aligner()
        #load geometry into thematic dictionary
        shapely_geom = shape(geometry)
        id = 1
        thematic_dict = {id: shapely_geom}
        loader = DictLoader(
            data_dict=thematic_dict
        )
        aligner.load_thematic_data(loader)
        #load reference data
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        dict_evaluated, prop_dictionary = aligner.evaluate(relevant_distances=relevant_distances)
        fc = aligner.get_results_as_geojson(
            resulttype=AlignerResultType.EVALUATED_PREDICTIONS,
            formula=True,
            attributes=False,
        )
        #return the featureclasses with the evaluated predictions
        return fc
    except Exception as e:
        return Response(str(e), status=400)

if __name__ == "__main__":
    with Configurator() as config:
        config.add_route("home", "/")
        config.add_route("actualiser", "/actualiser")
        config.scan()
        app = config.make_wsgi_app()
    server = make_server(host, port, app)
    print("webservice active on http://localhost:" + str(port))
    server.serve_forever()
