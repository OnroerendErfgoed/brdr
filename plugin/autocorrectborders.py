# -*- coding: utf-8 -*-

"""
***************************************************************************
*   version: v0.9.3 (proof of concept - working version)
*   author: Karel Dieussaert
*   history:
*            -initial version based on pyQGIS
*            -added exclusion of circles
*            -more efficient merge/union-logical
*            -removed resulting group layer (to prevent crashing of QGIS) -
                    extra research needed
*            -add logic for openbaar domein (od_strategy)
*            -intermediate layers added as an advanced parameter
*            -Native processes as child_algorithms
*            -Process NonThreaded to fix QGIS from crashing
*            -Added advanced parameter for processing input-multipolygons as single
                    polygons
*            -rewriting to use AutoReferencer (shapely-python)
*            -cleanup and added docs to AutoReferencer
*            -resulting output made available for further QGIS-modelling
*            -added enum - parameter to download actual GRB (adp-gbg-knw)
*            -added enum - parameter for od-strategy
*
MIT LICENSE:
Copyright (c) 2023-2024 Flanders Heritage Agency

Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify,
merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the
following conditions:

The above copyright notice and this permission notice shall be included in all copies
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***************************************************************************"""

from qgis import processing
from qgis.PyQt.QtCore import QCoreApplication
from qgis.PyQt.QtCore import QVariant
from qgis.PyQt.QtCore import Qt
from qgis.core import QgsCoordinateReferenceSystem
from qgis.core import QgsFeature
from qgis.core import QgsFeatureSink
from qgis.core import QgsField
from qgis.core import QgsGeometry
from qgis.core import QgsProcessing
from qgis.core import QgsProcessingAlgorithm
from qgis.core import QgsProcessingException
from qgis.core import QgsProcessingMultiStepFeedback
from qgis.core import QgsProcessingOutputVectorLayer
from qgis.core import QgsProcessingParameterBoolean
from qgis.core import QgsProcessingParameterEnum
from qgis.core import QgsProcessingParameterFeatureSource
from qgis.core import QgsProcessingParameterField
from qgis.core import QgsProcessingParameterNumber
from qgis.core import QgsProject
from qgis.core import QgsStyle
from qgis.core import QgsVectorLayer
from shapely import Polygon
from shapely import from_wkt
from shapely import to_wkt
from shapely import unary_union

try:
    import brdr

    if brdr.__version__ != "0.1.0":
        raise ValueError("Version mismatch")

except (ModuleNotFoundError, ValueError):
    print("Module package_name not found. Installing from PyPi")
    from pip._internal import main as pip

    pip(["install", "brdr==0.1.0"])
    import brdr

    print(brdr.__version__)

from brdr.auto_referencer import AutoReferencer


class AutocorrectBordersProcessingAlgorithm(QgsProcessingAlgorithm):
    """
    This script searches for overlap relevance between thematic borders and reference
    borders, and creates a resulting border based on the overlapping areas that are
    relevant.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT_THEMATIC = "INPUT_THEMATIC"
    INPUT_REFERENCE = "INPUT_REFERENCE"
    ENUM_REFERENCE = "ENUM_REFERENCE"
    ENUM_REFERENCE_OPTIONS = [
        "LOCAL REFERENCE LAYER (choose LAYER and ID below)",
        "download GRB - actuele percelen (adp)",
        "download GRB - actuele gebouwen (gbg)",
        "download GRB - actuele kunstwerken (knw)",
    ]
    SELECTED_REFERENCE = "None"
    ENUM_OD_STRATEGY = 'ENUM_OD_STRATEGY'
    ENUM_OD_STRATEGY_OPTIONS = [
        'EXCLUDE',
        'AS IS',
        'SNAP - ONE SIDE',
        'SNAP - ALL SIDE',
        'SNAP - BIG AREA'
    ]
    SELECTED_OD_STRATEGY = 'SNAP - ONE SIDE'
    RESULT = "RESULT"
    RESULT_DIFF = "RESULT_DIFF"
    RESULT_DIFF_PLUS = "RESULT_DIFF_PLUS"
    RESULT_DIFF_MIN = "RESULT_DIFF_MIN"
    OUTPUT_RESULT = "OUTPUT_RESULT"
    OUTPUT_RESULT_DIFF = "OUTPUT_RESULT_DIFF"
    OUTPUT_RESULT_DIFF_PLUS = "OUTPUT_RESULT_DIFF_PLUS"
    OUTPUT_RESULT_DIFF_MIN = "OUTPUT_RESULT_DIFF_MIN"

    INTERMEDIATE_LAYER_GROUP = "INTERMEDIATE_LAYER_GROUP"

    LAYER_RESULT = "LAYER_RESULT"
    LAYER_RESULT_DIFF = "LAYER_RESULT_DIFF"
    LAYER_RESULT_DIFF_PLUS = "LAYER_RESULT_DIFF_PLUS"
    LAYER_RESULT_DIFF_MIN = "LAYER_RESULT_DIFF_MIN"
    LAYER_SIGNIFICANT_INTERSECTION = "LAYER_SIGNIFICANT_INTERSECTION"
    LAYER_SIGNIFICANT_DIFFERENCE = "LAYER_SIGNIFICANT_DIFFERENCE"
    LAYER_REFERENCE = "LAYER_REFERENCE"

    SUFFIX = ""
    # theme_ID (can be a multipolygon)
    ID_THEME_GLOBAL = "id_theme"
    # theme_ID for singlified multipolygons following syntax #theme_id + '_' + row_nr
    ID_THEME = "id_theme_singlified"
    ID_REFERENCE = "id_ref"
    OVERLAY_FIELDS_PREFIX = ""
    OD_STRATEGY = 0
    THRESHOLD_OVERLAP_PERCENTAGE = 50
    THRESHOLD_EXCLUSION_AREA = 0
    THRESHOLD_EXCLUSION_PERCENTAGE = 0
    RELEVANT_DISTANCE = 0
    BUFFER_DISTANCE = 0
    THRESHOLD_CIRCLE_RATIO = 0.98
    CORR_DISTANCE = 0.01
    SHOW_INTERMEDIATE_LAYERS = False
    PROCESS_MULTI_AS_SINGLE_POLYGONS = True
    MITRE_LIMIT = 10
    CRS = "EPSG:31370"
    QUAD_SEGS = 5
    BUFFER_MULTIPLICATION_FACTOR = 1.01
    DOWNLOAD_LIMIT = 10000
    MAX_REFERENCE_BUFFER = 10

    def flags(self):
        return super().flags() | QgsProcessingAlgorithm.FlagNoThreading

    @staticmethod
    def tr(string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate("Processing", string)

    def createInstance(self):
        return AutocorrectBordersProcessingAlgorithm()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return "autocorrectBorders"

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr("Autocorrectborders")

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr("autocorrectborders - scripts")

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return "autocorrectbordersscripts"

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it.
        """
        return self.tr(
            "This script searches for overlap relevance between thematic borders and "
            "reference borders, and creates a resulting border based on the overlapping "
            "areas that are relevant"
        )

    def geom_shapely_to_qgis(self, geom_shapely):
        """
        Method to convert a Shapely-geometry to a QGIS geometry
        """
        wkt = to_wkt(geom_shapely, rounding_precision=-1, output_dimension=2)
        geom_qgis = QgsGeometry.fromWkt(wkt)
        return geom_qgis

    def geom_qgis_to_shapely(self, geom_qgis):
        """
        Method to convert a QGIS-geometry to a Shapely-geometry
        """
        wkt = geom_qgis.asWkt()
        geom_shapely = from_wkt(wkt)
        return geom_shapely

    def create_temp_layer(
        self, name, group_layer_name, field_name, style_name, visible
    ):
        """
        Create a temporary QGIS layer in the TOC based on:
         name: name for the temporary layer
         field_name: the ID-fieldname that has to be added
         stylename: name of a predefined layer-style
         visible: if the layer has to be checked or unchecked
        """
        # https://docs.qgis.org/3.28/en/docs/pyqgis_developer_cookbook/cheat_sheet.html
        qinst = QgsProject.instance()
        lyrs = qinst.mapLayersByName(name)
        root = qinst.layerTreeRoot()
        if len(lyrs) != 0:
            for lyr in lyrs:
                root.removeLayer(lyr)
                qinst.removeMapLayer(lyr.id())
        vl = QgsVectorLayer("MultiPolygon", name, "memory")
        vl.setCrs(QgsCoordinateReferenceSystem(self.CRS))
        pr = vl.dataProvider()
        pr.addAttributes([QgsField(field_name, QVariant.String)])
        vl.updateFields()
        # styling
        vl.setOpacity(0.5)
        if style_name is None or style_name == "":
            symbol = None
        else:
            symbol = QgsStyle.defaultStyle().symbol(style_name)
        if symbol is not None:
            vl.renderer().setSymbol(symbol)
        # adding layer to TOC
        qinst.addMapLayer(
            vl, False
        )  # False so that it doesn't get inserted at default position
        root.insertLayer(0, vl)
        node = root.findLayer(vl.id())
        if node:
            new_state = Qt.Checked if visible else Qt.Unchecked
            node.setItemVisibilityChecked(new_state)
        vl.triggerRepaint()
        return vl

    def get_layer_by_name(self, layer_name):
        """
        Get the layer-object based on the layername
        """
        layers = QgsProject.instance().mapLayersByName(layer_name)
        return layers[0]

    def add_geom_to_temp_layer(self, layer_name, geom, id):
        """
        Method to add a feature to a temporary layer
        layer_name: name of the temporary layer where the feature has to be added
        geom: geometry that has to be added
        id: ID related to the geometry
        """
        layer = self.get_layer_by_name(layer_name)
        f = QgsFeature()
        f.setGeometry(geom)
        f.setAttributes([id])
        pr = layer.dataProvider()
        pr.addFeature(f)
        layer.updateExtents()
        layer.triggerRepaint()
        return

    def merge_multi_to_single(self, dictionary):
        """
        Merges geometries in a dictionary from multiple themes into a single theme.

        Args: dictionary (dict): A dictionary where keys are theme IDs and values are
            geometries (e.g., shapely Polygon objects).

        Returns: dict: A new dictionary with merged geometries, where keys are global
            theme IDs and values are merged geometries.

        """
        dict_out = {}
        for id_theme in dictionary:
            id_theme_global = id_theme.split("_")[0]
            geom = dictionary[id_theme]
            if geom.is_empty or geom is None:
                continue
            arr = [geom]
            if id_theme_global not in dict_out:
                dict_out[id_theme_global] = [Polygon()]
            lst = dict_out[id_theme_global]
            lst.extend(arr)
            dict_out[id_theme_global] = lst
        for id_theme_global in dict_out:
            dict_out[id_theme_global] = unary_union(dict_out[id_theme_global])
        return dict_out

    def geom_from_dict(self, dict, key):
        """
        Get the geometry from a dictionary with geometries. If key not present,
        an empty Polygon is returned
        """
        if key in dict:
            geom = dict[key]
        else:
            geom = Polygon()
        return geom

    def layer_to_featuresink(self, parameters, context, layer_name, sink_name):
        layer = self.get_layer_by_name(layer_name)
        source = layer.dataProvider()
        (sink, dest_id) = self.parameterAsSink(
            parameters,
            sink_name,
            context,
            source.fields(),
            source.wkbType(),
            QgsCoordinateReferenceSystem(self.CRS),
        )
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))
        features = source.getFeatures()
        for current, feature in enumerate(features):
            sink.addFeature(feature, QgsFeatureSink.FastInsert)
        return sink, dest_id

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # standard parameters
        parameter = QgsProcessingParameterFeatureSource(
            self.INPUT_THEMATIC,
            self.tr("THEMATIC LAYER"),
            [QgsProcessing.TypeVectorAnyGeometry],
            defaultValue="themelayer",
        )
        parameter.setFlags(parameter.flags())
        self.addParameter(parameter)
        parameter = QgsProcessingParameterField(
            self.ID_THEME_GLOBAL,
            "Choose thematic ID",
            "theme_identifier",
            self.INPUT_THEMATIC,
        )
        parameter.setFlags(parameter.flags())
        self.addParameter(parameter)
        parameter = QgsProcessingParameterEnum(
            self.ENUM_REFERENCE,
            "Select Reference Layer:",
            options=self.ENUM_REFERENCE_OPTIONS,
            defaultValue=0,  # Index of the default option (e.g., 'Option A')
        )
        parameter.setFlags(parameter.flags())
        self.addParameter(parameter)
        parameter = QgsProcessingParameterFeatureSource(
            self.INPUT_REFERENCE,
            self.tr("REFERENCE LAYER"),
            [QgsProcessing.TypeVectorAnyGeometry],
            defaultValue="referencelayer",
        )
        parameter.setFlags(parameter.flags())
        self.addParameter(parameter)
        parameter = QgsProcessingParameterField(
            self.ID_REFERENCE,
            "Choose reference ID",
            "ref_identifier",
            self.INPUT_REFERENCE,
        )
        parameter.setFlags(parameter.flags())
        self.addParameter(parameter)
        parameter = QgsProcessingParameterNumber(
            "RELEVANT_DISTANCE",
            "RELEVANT_DISTANCE (meter)",
            type=QgsProcessingParameterNumber.Double,
            defaultValue=2,
        )
        parameter.setFlags(parameter.flags())
        self.addParameter(parameter)
        parameter = QgsProcessingParameterEnum(
            self.ENUM_OD_STRATEGY,
            'Select OD-STRATEGY:',
            options=self.ENUM_OD_STRATEGY_OPTIONS,
            defaultValue=2  # Index of the default option (e.g., 'Snap - one side')
        )
        parameter.setFlags(parameter.flags())
        self.addParameter(parameter)

        self.addOutput(
            QgsProcessingOutputVectorLayer(
                self.OUTPUT_RESULT,
                self.LAYER_RESULT,
                QgsProcessing.TypeVectorAnyGeometry,
            )
        )
        self.addOutput(
            QgsProcessingOutputVectorLayer(
                self.OUTPUT_RESULT_DIFF,
                self.LAYER_RESULT_DIFF,
                QgsProcessing.TypeVectorAnyGeometry,
            )
        )
        self.addOutput(
            QgsProcessingOutputVectorLayer(
                self.OUTPUT_RESULT_DIFF_PLUS,
                self.LAYER_RESULT_DIFF_PLUS,
                QgsProcessing.TypeVectorAnyGeometry,
            )
        )
        self.addOutput(
            QgsProcessingOutputVectorLayer(
                self.OUTPUT_RESULT_DIFF_MIN,
                self.LAYER_RESULT_DIFF_MIN,
                QgsProcessing.TypeVectorAnyGeometry,
            )
        )
        # advanced parameters
        parameter = QgsProcessingParameterNumber(
            "THRESHOLD_OVERLAP_PERCENTAGE",
            "THRESHOLD_OVERLAP_PERCENTAGE (%)",
            type=QgsProcessingParameterNumber.Double,
            defaultValue=50,
        )
        self.addParameter(parameter)
        parameter = QgsProcessingParameterBoolean(
            "PROCESS_MULTI_AS_SINGLE_POLYGONS",
            "PROCESS_MULTI_AS_SINGLE_POLYGONS",
            defaultValue=True,
        )
        # parameter.setFlags(parameter.flags() |
        # QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(parameter)
        parameter = QgsProcessingParameterBoolean(
            "SHOW_INTERMEDIATE_LAYERS", "SHOW_INTERMEDIATE_LAYERS", defaultValue=False
        )
        # parameter.setFlags(parameter.flags() |
        # QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(parameter)

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """

        feedback_steps = 6
        feedback = QgsProcessingMultiStepFeedback(feedback_steps, feedback)
        feedback.pushInfo("START")
        outputs = {}

        self.prepare_parameters(parameters)

        thematic, thematic_buffered = self._thematic_preparation(
            context, feedback, outputs, parameters
        )
        if thematic is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.test))

        # Load thematic into a shapely_dict:
        dict_thematic = {}
        features = thematic.getFeatures()
        for current, feature in enumerate(features):
            if feedback.isCanceled():
                return {}
            id_theme = feature.attribute(self.ID_THEME)
            dict_thematic[id_theme] = self.geom_qgis_to_shapely(feature.geometry())
        feedback.pushInfo("1) BEREKENING - Thematic layer fixed")
        feedback.setCurrentStep(1)
        if feedback.isCanceled():
            return {}

        # REFERENCE PREPARATION
        if self.SELECTED_REFERENCE == 0:
            outputs[self.INPUT_REFERENCE + "_extract"] = processing.run(
                "native:extractbylocation",
                {
                    "INPUT": parameters[self.INPUT_REFERENCE],
                    "PREDICATE": [0],
                    "INTERSECT": thematic_buffered,
                    "OUTPUT": "TEMPORARY_OUTPUT",
                },
                context=context,
                feedback=feedback,
                is_child_algorithm=True,
            )
            reference = context.getMapLayer(
                outputs[self.INPUT_REFERENCE + "_extract"]["OUTPUT"]
            )
            if reference.sourceCrs().authid() != self.CRS:
                raise QgsProcessingException(
                    "Thematic layer and ReferenceLayer are in a different CRS. "
                    "Please provide them in the same CRS (EPSG:31370 or EPSG:3812)"
                )
            outputs[self.INPUT_REFERENCE + "_id"] = processing.run(
                "native:fieldcalculator",
                {
                    "INPUT": reference,
                    "FIELD_NAME": self.ID_REFERENCE,
                    "FIELD_TYPE": 2,
                    "FIELD_LENGTH": 0,
                    "FIELD_PRECISION": 0,
                    "FORMULA": "to_string(" + parameters[self.ID_REFERENCE] + ")",
                    "OUTPUT": "TEMPORARY_OUTPUT",
                },
                context=context,
                feedback=feedback,
                is_child_algorithm=True,
            )
            reference = context.getMapLayer(
                outputs[self.INPUT_REFERENCE + "_id"]["OUTPUT"]
            )
            outputs[self.INPUT_REFERENCE + "_fixed"] = processing.run(
                "native:fixgeometries",
                {
                    "INPUT": reference,
                    "METHOD": 1,
                    "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT,
                },
                context=context,
                feedback=feedback,
                is_child_algorithm=True,
            )
            reference = context.getMapLayer(
                outputs[self.INPUT_REFERENCE + "_fixed"]["OUTPUT"]
            )
            outputs[self.INPUT_REFERENCE + "_dropMZ"] = processing.run(
                "native:dropmzvalues",
                {
                    "INPUT": reference,
                    "DROP_M_VALUES": True,
                    "DROP_Z_VALUES": True,
                    "OUTPUT": "TEMPORARY_OUTPUT",
                },
                context=context,
                feedback=feedback,
                is_child_algorithm=True,
            )
            reference = context.getMapLayer(
                outputs[self.INPUT_REFERENCE + "_dropMZ"]["OUTPUT"]
            )
            if reference is None:
                raise QgsProcessingException(
                    self.invalidSourceError(parameters, self.INPUT_REFERENCE)
                )

            # Load reference into a shapely_dict:
            dict_reference = {}
            features = reference.getFeatures()
            for current, feature in enumerate(features):
                if feedback.isCanceled():
                    return {}
                id_reference = feature.attribute(self.ID_REFERENCE)
                dict_reference[id_reference] = self.geom_qgis_to_shapely(
                    feature.geometry()
                )
            feedback.pushInfo("2) BEREKENING - Reference layer fixed")
            feedback.setCurrentStep(2)
            if feedback.isCanceled():
                return {}

        # MAKE TEMPORARY LAYERS
        if self.SELECTED_REFERENCE != 0:
            self.create_temp_layer(
                self.LAYER_REFERENCE,
                self.INTERMEDIATE_LAYER_GROUP,
                self.ID_REFERENCE,
                "gray 1 fill",
                True,
            )
        if self.SHOW_INTERMEDIATE_LAYERS:
            self.create_temp_layer(
                self.LAYER_SIGNIFICANT_INTERSECTION,
                self.INTERMEDIATE_LAYER_GROUP,
                self.ID_THEME_GLOBAL,
                "simple green fill",
                False,
            )
            self.create_temp_layer(
                self.LAYER_SIGNIFICANT_DIFFERENCE,
                self.INTERMEDIATE_LAYER_GROUP,
                self.ID_THEME_GLOBAL,
                "simple red fill",
                False,
            )
        self.create_temp_layer(
            self.LAYER_RESULT,
            self.INTERMEDIATE_LAYER_GROUP,
            self.ID_THEME_GLOBAL,
            "outline xpattern",
            True,
        )
        self.create_temp_layer(
            self.LAYER_RESULT_DIFF,
            self.INTERMEDIATE_LAYER_GROUP,
            self.ID_THEME_GLOBAL,
            "hashed black X",
            False,
        )
        self.create_temp_layer(
            self.LAYER_RESULT_DIFF_PLUS,
            self.INTERMEDIATE_LAYER_GROUP,
            self.ID_THEME_GLOBAL,
            "hashed cgreen /",
            False,
        )
        self.create_temp_layer(
            self.LAYER_RESULT_DIFF_MIN,
            self.INTERMEDIATE_LAYER_GROUP,
            self.ID_THEME_GLOBAL,
            "hashed cred /",
            False,
        )

        # AUTOREFERENCER IMPLEMENTATION
        auto_referencer = AutoReferencer(
            feedback=feedback,
            relevant_distance=self.RELEVANT_DISTANCE,
            threshold_overlap_percentage=self.THRESHOLD_OVERLAP_PERCENTAGE,
        )

        # set parameters
        auto_referencer.feedback = feedback
        auto_referencer.relevant_distance = self.RELEVANT_DISTANCE
        auto_referencer.od_strategy = self.OD_STRATEGY
        auto_referencer.THRESHOLD_CIRCLE_RATIO = self.THRESHOLD_CIRCLE_RATIO
        auto_referencer.THRESHOLD_EXCLUSION_AREA = self.THRESHOLD_EXCLUSION_AREA
        auto_referencer.THRESHOLD_EXCLUSION_PERCENTAGE = (
            self.THRESHOLD_EXCLUSION_PERCENTAGE
        )
        auto_referencer.CORR_DISTANCE = self.CORR_DISTANCE
        auto_referencer.MITRE_LIMIT = self.MITRE_LIMIT
        auto_referencer.QUAD_SEGS = self.QUAD_SEGS
        auto_referencer.BUFFER_MULTIPLICATION_FACTOR = self.BUFFER_MULTIPLICATION_FACTOR
        auto_referencer.MAX_REFERENCE_BUFFER = self.MAX_REFERENCE_BUFFER
        auto_referencer.CRS = self.CRS
        auto_referencer.DOWNLOAD_LIMIT = self.DOWNLOAD_LIMIT

        feedback.pushInfo("Load thematic data")
        auto_referencer.load_thematic_data_dict(dict_thematic)

        feedback.pushInfo("Load reference data")
        if self.SELECTED_REFERENCE == 0:
            auto_referencer.load_reference_data_dict(dict_reference)
        elif self.SELECTED_REFERENCE == 1:  # adp
            auto_referencer.load_reference_data_grb_actual()
        elif self.SELECTED_REFERENCE == 2:  # gbg
            auto_referencer.load_reference_data_grb_actual(grb_type="gbg")
        elif self.SELECTED_REFERENCE == 3:  # knw
            auto_referencer.load_reference_data_grb_actual(grb_type="knw")

        #
        feedback.pushInfo("START PROCESSING")
        (
            dict_result,
            dict_result_diff,
            dict_result_diff_plus,
            dict_result_diff_min,
            dict_relevant_intersection,
            dict_relevant_diff,
        ) = auto_referencer.process_dict_thematic(
            self.RELEVANT_DISTANCE, self.OD_STRATEGY, self.THRESHOLD_OVERLAP_PERCENTAGE
        )

        feedback.pushInfo("END PROCESSING")
        # write results to output-layers
        feedback.pushInfo("Merging multi_to_single RESULTS")
        result_merged = self.merge_multi_to_single(dict_result)
        result_diff_merged = self.merge_multi_to_single(dict_result_diff)
        result_diff_plus_plus_merged = self.merge_multi_to_single(dict_result_diff_plus)
        result_diff_min_min_merged = self.merge_multi_to_single(dict_result_diff_min)
        relevant_intersection_merged = self.merge_multi_to_single(
            dict_relevant_intersection
        )
        relevant_diff_merged = self.merge_multi_to_single(dict_relevant_diff)

        # write results to output-layers
        feedback.pushInfo("WRITING RESULTS")
        if self.SELECTED_REFERENCE != 0:
            dict_reference = auto_referencer.dict_reference
            for id_reference in dict_reference:
                self.add_geom_to_temp_layer(
                    self.LAYER_REFERENCE,
                    self.geom_shapely_to_qgis(
                        self.geom_from_dict(dict_reference, id_reference)
                    ),
                    id_reference,
                )

        list_id_theme_global = []
        for id_theme in dict_thematic:
            id_theme_global = id_theme.split("_")[0]
            if id_theme_global not in list_id_theme_global:
                list_id_theme_global.append(id_theme_global)
        for id_theme_global in list_id_theme_global:
            self.add_geom_to_temp_layer(
                self.LAYER_RESULT,
                self.geom_shapely_to_qgis(
                    self.geom_from_dict(result_merged, id_theme_global)
                ),
                id_theme_global,
            )
            self.add_geom_to_temp_layer(
                self.LAYER_RESULT_DIFF,
                self.geom_shapely_to_qgis(
                    self.geom_from_dict(result_diff_merged, id_theme_global)
                ),
                id_theme_global,
            )
            self.add_geom_to_temp_layer(
                self.LAYER_RESULT_DIFF_PLUS,
                self.geom_shapely_to_qgis(
                    self.geom_from_dict(result_diff_plus_plus_merged, id_theme_global)
                ),
                id_theme_global,
            )
            self.add_geom_to_temp_layer(
                self.LAYER_RESULT_DIFF_MIN,
                self.geom_shapely_to_qgis(
                    self.geom_from_dict(result_diff_min_min_merged, id_theme_global)
                ),
                id_theme_global,
            )
            if self.SHOW_INTERMEDIATE_LAYERS:
                self.add_geom_to_temp_layer(
                    self.LAYER_SIGNIFICANT_INTERSECTION,
                    self.geom_shapely_to_qgis(
                        self.geom_from_dict(
                            relevant_intersection_merged, id_theme_global
                        )
                    ),
                    id_theme_global,
                )
                self.add_geom_to_temp_layer(
                    self.LAYER_SIGNIFICANT_DIFFERENCE,
                    self.geom_shapely_to_qgis(
                        self.geom_from_dict(relevant_diff_merged, id_theme_global)
                    ),
                    id_theme_global,
                )
        self.RESULT = QgsProject.instance().mapLayersByName(self.LAYER_RESULT)[0]
        self.RESULT_DIFF = QgsProject.instance().mapLayersByName(
            self.LAYER_RESULT_DIFF
        )[0]
        self.RESULT_DIFF_PLUS = QgsProject.instance().mapLayersByName(
            self.LAYER_RESULT_DIFF_PLUS
        )[0]
        self.RESULT_DIFF_MIN = QgsProject.instance().mapLayersByName(
            self.LAYER_RESULT_DIFF_MIN
        )[0]

        QgsProject.instance().reloadAllLayers()

        feedback.pushInfo("Resulterende geometrie berekend")
        feedback.setCurrentStep(6)
        if feedback.isCanceled():
            return {}

        feedback.pushInfo("EINDE: RESULTAAT BEREKEND")
        return {
            self.OUTPUT_RESULT: self.RESULT,
            self.OUTPUT_RESULT_DIFF: self.RESULT_DIFF,
            self.OUTPUT_RESULT_DIFF_PLUS: self.RESULT_DIFF_PLUS,
            self.OUTPUT_RESULT_DIFF_MIN: self.RESULT_DIFF_MIN,
        }

    def _thematic_preparation(self, context, feedback, outputs, parameters):
        # THEMATIC PREPARATION
        outputs[self.INPUT_THEMATIC + "_id"] = processing.run(
            "native:fieldcalculator",
            {
                "INPUT": parameters[self.INPUT_THEMATIC],
                "FIELD_NAME": self.ID_THEME_GLOBAL,
                "FIELD_TYPE": 2,
                "FIELD_LENGTH": 0,
                "FIELD_PRECISION": 0,
                "FORMULA": "to_string(" + parameters[self.ID_THEME_GLOBAL] + ")",
                "OUTPUT": "TEMPORARY_OUTPUT",
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )
        thematic = context.getMapLayer(outputs[self.INPUT_THEMATIC + "_id"]["OUTPUT"])
        self.CRS = (
            thematic.sourceCrs().authid()
        )  # set CRS for the calculations, based on the THEMATIC input layer
        if self.PROCESS_MULTI_AS_SINGLE_POLYGONS:
            outputs[self.INPUT_THEMATIC + "_single"] = processing.run(
                "native:multiparttosingleparts",
                {"INPUT": thematic, "OUTPUT": "TEMPORARY_OUTPUT"},
                context=context,
                feedback=feedback,
                is_child_algorithm=True,
            )
            thematic = context.getMapLayer(
                outputs[self.INPUT_THEMATIC + "_single"]["OUTPUT"]
            )
        outputs[self.INPUT_THEMATIC + "_single_id"] = processing.run(
            "native:fieldcalculator",
            {
                "INPUT": thematic,
                "FIELD_NAME": self.ID_THEME,
                "FIELD_TYPE": 2,
                "FIELD_LENGTH": 0,
                "FIELD_PRECISION": 0,
                "FORMULA": "to_string("
                + parameters[self.ID_THEME_GLOBAL]
                + ") + '_'+ to_string(@id)",
                "OUTPUT": "TEMPORARY_OUTPUT",
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )
        thematic = context.getMapLayer(
            outputs[self.INPUT_THEMATIC + "_single_id"]["OUTPUT"]
        )
        outputs[self.INPUT_THEMATIC + "_fixed"] = processing.run(
            "native:fixgeometries",
            {"INPUT": thematic, "METHOD": 1, "OUTPUT": QgsProcessing.TEMPORARY_OUTPUT},
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )
        thematic = context.getMapLayer(
            outputs[self.INPUT_THEMATIC + "_fixed"]["OUTPUT"]
        )
        outputs[self.INPUT_THEMATIC + "_enriched"] = processing.run(
            "qgis:exportaddgeometrycolumns",
            {"INPUT": thematic, "CALC_METHOD": 0, "OUTPUT": "TEMPORARY_OUTPUT"},
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )
        thematic = context.getMapLayer(
            outputs[self.INPUT_THEMATIC + "_enriched"]["OUTPUT"]
        )
        outputs[self.INPUT_THEMATIC + "_dropMZ"] = processing.run(
            "native:dropmzvalues",
            {
                "INPUT": thematic,
                "DROP_M_VALUES": True,
                "DROP_Z_VALUES": True,
                "OUTPUT": "TEMPORARY_OUTPUT",
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )
        thematic = context.getMapLayer(
            outputs[self.INPUT_THEMATIC + "_dropMZ"]["OUTPUT"]
        )
        # buffer the thematic layer to select all plots around it that are relevant to
        # the calculations
        outputs[self.INPUT_THEMATIC + "_buffered"] = processing.run(
            "native:buffer",
            {
                "INPUT": thematic,
                "DISTANCE": self.BUFFER_MULTIPLICATION_FACTOR * self.RELEVANT_DISTANCE,
                "SEGMENTS": self.QUAD_SEGS,
                "END_CAP_STYLE": 0,
                "JOIN_STYLE": 1,
                "MITRE_LIMIT": self.MITRE_LIMIT,
                "DISSOLVE": False,
                "OUTPUT": "TEMPORARY_OUTPUT",
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )
        thematic_buffered = context.getMapLayer(
            outputs[self.INPUT_THEMATIC + "_buffered"]["OUTPUT"]
        )
        return thematic, thematic_buffered

    def prepare_parameters(self, parameters):
        # PARAMETER PREPARATION
        self.RELEVANT_DISTANCE = parameters["RELEVANT_DISTANCE"]
        self.BUFFER_DISTANCE = self.RELEVANT_DISTANCE / 2
        self.THRESHOLD_OVERLAP_PERCENTAGE = parameters["THRESHOLD_OVERLAP_PERCENTAGE"]
        self.SELECTED_OD_STRATEGY = parameters[self.ENUM_OD_STRATEGY]
        self.OD_STRATEGY = self.SELECTED_OD_STRATEGY - 1
        self.SHOW_INTERMEDIATE_LAYERS = parameters["SHOW_INTERMEDIATE_LAYERS"]
        self.PROCESS_MULTI_AS_SINGLE_POLYGONS = parameters[
            "PROCESS_MULTI_AS_SINGLE_POLYGONS"
        ]
        self.SUFFIX = "_" + str(self.RELEVANT_DISTANCE) + "_OD_" + str(self.OD_STRATEGY)
        self.LAYER_SIGNIFICANT_INTERSECTION = (
            self.LAYER_SIGNIFICANT_INTERSECTION + self.SUFFIX
        )
        self.LAYER_SIGNIFICANT_DIFFERENCE = (
            self.LAYER_SIGNIFICANT_DIFFERENCE + self.SUFFIX
        )
        self.LAYER_RESULT = self.LAYER_RESULT + self.SUFFIX
        self.LAYER_RESULT_DIFF = self.LAYER_RESULT_DIFF + self.SUFFIX
        self.LAYER_RESULT_DIFF_PLUS = self.LAYER_RESULT_DIFF_PLUS + self.SUFFIX
        self.LAYER_RESULT_DIFF_MIN = self.LAYER_RESULT_DIFF_MIN + self.SUFFIX
        self.SELECTED_REFERENCE = parameters[self.ENUM_REFERENCE]
        self.LAYER_REFERENCE = self.ENUM_REFERENCE_OPTIONS[
            parameters[self.ENUM_REFERENCE]
        ]
