from enum import Enum

import requests
from requests import RequestException


class BRKCollectionLoader:
    @classmethod
    def _fetch_values(cls):
        _url = "https://api.pdok.nl/kadaster/brk-kadastrale-kaart/ogc/v1/collections?f=json"
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
                "kadastralegrens": "KadastraleGrens",
                "perceel": "Perceel",
                "openbareruimtenaam": "OpenbareRuimteNaam",
                "bebouwing": "Bebouwing",
                "nummeraanduidingreeks": "Nummeraanduidingreeks",
            }
        return Enum("BRKType", dict_values)


BRKType = BRKCollectionLoader.get_enum()
