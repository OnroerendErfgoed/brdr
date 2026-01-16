from enum import Enum


class OEType(str, Enum):
    """
    Different types of Onroerend Eefgoed-objects are available:

    * AO: aanduidingsobjecten
    * EO: erfgoedobjecten
    """

    AO = "aanduidingsobjecten"
    EO = "erfgoedobjecten"
