from enum import Enum

import requests
from requests import RequestException


class GRBTypeLoader:
    @classmethod
    def _fetch_values(cls):
        _url = (
            "https://geo.api.vlaanderen.be/GRB/ogc/features/collections"
            + "/?f=application%2Fjson"
        )
        response = requests.get(_url, timeout=30)
        response.raise_for_status()
        data = response.json()
        dict_values = {}
        for coll in data["collections"]:
            dict_values[coll["id"]] = coll["title"]
        return dict_values

    @classmethod
    def get_enum(cls):
        try:
            dict_values = cls._fetch_values()
        except (RequestException, ValueError, KeyError):
            dict_values = {
                "ADP": "Administratieve percelen",
                "GBG": "gebouwen",
                "KNW": "kunstwerken",
                "WBN": "wegbaan",
                "WGO": "wegopdeling",
            }
        return Enum("GRBType", dict_values)


GRBType = GRBTypeLoader.get_enum()
