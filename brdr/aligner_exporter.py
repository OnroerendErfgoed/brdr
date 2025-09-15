"""
Aligner Exporter - handles export and serialization functionality
"""
import os
from typing import Dict, Any, Optional

from brdr.aligner_core import AlignerCore
from brdr.enums import AlignerResultType, AlignerInputType
from brdr.typings import ProcessResult
from brdr.utils import (
    get_dict_geojsons_from_series_dict,
    write_geojson,
    geojson_from_dict
)


class AlignerExporter:
    """
    Exporter class responsible for serialization and export functionality.
    
    This class handles:
    - Exporting results as GeoJSON
    - Saving results to files
    - Formatting output data
    - Managing export configurations
    """
    
    def __init__(self, core: AlignerCore):
        """
        Initialize the exporter with a core aligner instance.
        
        Args:
            core: AlignerCore instance with loaded data
        """
        self.core = core
        self.config = core.config
        self.logger = core.logger
    
    def get_results_as_geojson(
        self,
        dict_processresults: Dict[Any, Dict[float, ProcessResult]] = None,
        resulttype: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        formula: bool = False,
        attributes: bool = True,
    ) -> Dict[str, Any]:
        """
        Export processing results as GeoJSON feature collections.
        
        Args:
            dict_processresults: Processing results to export (None for all)
            resulttype: Type of results to export
            formula: Whether to include formula information
            attributes: Whether to include geometry attributes
            
        Returns:
            Dictionary of GeoJSON feature collections by result type
            
        Raises:
            ValueError: If no results available or invalid parameters
        """
        if dict_processresults is None:
            raise ValueError("No process results provided")
        
        if not dict_processresults:
            self.logger.feedback_warning("No process results to export")
            return {}
        
        self.logger.feedback_info(f"Exporting results as GeoJSON (type: {resulttype})")
        
        # Get properties if needed
        series_prop_dict = None
        if formula or attributes:
            series_prop_dict = self._prepare_properties_for_export(
                dict_processresults, formula, attributes
            )
        
        # Convert to GeoJSON
        geojson_dict = get_dict_geojsons_from_series_dict(
            dict_processresults,
            self.config.crs,
            "brdr_id",
            series_prop_dict,
            geom_attributes=attributes
        )
        
        self.logger.feedback_info(f"Exported {len(geojson_dict)} GeoJSON collections")
        return geojson_dict
    
    def save_results(
        self,
        path: str,
        dict_processresults: Dict[Any, Dict[float, ProcessResult]] = None,
        resulttype: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        formula: bool = True,
        attributes: bool = True,
    ) -> None:
        """
        Save processing results to files.
        
        Args:
            path: Directory path to save files
            dict_processresults: Processing results to save
            resulttype: Type of results to save
            formula: Whether to include formula information
            attributes: Whether to include geometry attributes
            
        Raises:
            ValueError: If path is invalid or no results to save
            OSError: If unable to create directory or write files
        """
        if not path:
            raise ValueError("Path cannot be empty")
        
        if dict_processresults is None:
            raise ValueError("No process results provided")
        
        # Ensure directory exists
        try:
            os.makedirs(path, exist_ok=True)
        except OSError as e:
            raise OSError(f"Unable to create directory {path}: {e}")
        
        self.logger.feedback_info(f"Saving results to {path}")
        
        # Get GeoJSON results
        geojson_dict = self.get_results_as_geojson(
            dict_processresults, resulttype, formula, attributes
        )
        
        # Save each result type to a separate file
        saved_files = []
        for result_key, geojson_data in geojson_dict.items():
            filename = f"{resulttype.value}_{result_key}.geojson"
            filepath = os.path.join(path, filename)
            
            try:
                write_geojson(filepath, geojson_data)
                saved_files.append(filename)
                self.logger.feedback_debug(f"Saved {filename}")
            except Exception as e:
                self.logger.feedback_warning(f"Failed to save {filename}: {e}")
        
        self.logger.feedback_info(f"Saved {len(saved_files)} files to {path}")
    
    def get_input_as_geojson(
        self, 
        inputtype: AlignerInputType = AlignerInputType.REFERENCE,
        include_properties: bool = True
    ):
        """
        Get input data as GeoJSON FeatureCollection.
        
        Args:
            inputtype: Type of input data to export
            include_properties: Whether to include feature properties
            
        Returns:
            GeoJSON FeatureCollection of input data
        """
        return self.core.get_input_as_geojson(inputtype)
    
    def export_thematic_data(
        self, 
        filepath: str, 
        include_properties: bool = True,
        attributes: bool = True
    ) -> None:
        """
        Export thematic data to GeoJSON file.
        
        Args:
            filepath: Path to save the GeoJSON file
            include_properties: Whether to include feature properties
            attributes: Whether to include geometry attributes
            
        Raises:
            ValueError: If no thematic data loaded
            OSError: If unable to write file
        """
        if not self.core.dict_thematic:
            raise ValueError("No thematic data loaded")
        
        self.logger.feedback_info(f"Exporting thematic data to {filepath}")
        
        # Prepare properties
        properties_dict = None
        if include_properties:
            properties_dict = self.core.dict_thematic_properties
        
        # Create GeoJSON
        geojson = geojson_from_dict(
            self.core.dict_thematic,
            self.config.crs,
            "theme_id",
            properties_dict,
            geom_attributes=attributes
        )
        
        # Save to file
        try:
            write_geojson(filepath, geojson)
            self.logger.feedback_info(f"Thematic data exported successfully")
        except Exception as e:
            raise OSError(f"Failed to write file {filepath}: {e}")
    
    def export_reference_data(
        self, 
        filepath: str, 
        include_properties: bool = True,
        attributes: bool = True
    ) -> None:
        """
        Export reference data to GeoJSON file.
        
        Args:
            filepath: Path to save the GeoJSON file
            include_properties: Whether to include feature properties
            attributes: Whether to include geometry attributes
            
        Raises:
            ValueError: If no reference data loaded
            OSError: If unable to write file
        """
        if not self.core.dict_reference:
            raise ValueError("No reference data loaded")
        
        self.logger.feedback_info(f"Exporting reference data to {filepath}")
        
        # Prepare properties
        properties_dict = None
        if include_properties:
            properties_dict = self.core.dict_reference_properties
        
        # Create GeoJSON
        geojson = geojson_from_dict(
            self.core.dict_reference,
            self.config.crs,
            "ref_id",
            properties_dict,
            geom_attributes=attributes
        )
        
        # Save to file
        try:
            write_geojson(filepath, geojson)
            self.logger.feedback_info(f"Reference data exported successfully")
        except Exception as e:
            raise OSError(f"Failed to write file {filepath}: {e}")
    
    def get_diff_metrics(
        self,
        dict_processresults: Dict[Any, Dict[float, ProcessResult]] = None,
        dict_thematic: Dict = None,
    ) -> Dict:
        """
        Calculate difference metrics for processing results.
        
        Args:
            dict_processresults: Processing results to analyze
            dict_thematic: Original thematic data for comparison
            
        Returns:
            Dictionary containing difference metrics
        """
        if dict_processresults is None:
            raise ValueError("No process results provided")
        
        if dict_thematic is None:
            dict_thematic = self.core.dict_thematic
        
        self.logger.feedback_info("Calculating difference metrics")
        
        # TODO: Implement difference metrics calculation
        # This would analyze the differences between original and processed geometries
        
        metrics = {}
        for theme_id, results in dict_processresults.items():
            if theme_id not in dict_thematic:
                continue
            
            original_geom = dict_thematic[theme_id]
            theme_metrics = {}
            
            for distance, result in results.items():
                if 'result' not in result:
                    continue
                
                processed_geom = result['result']
                
                # Calculate basic metrics
                theme_metrics[distance] = self._calculate_single_metrics(
                    original_geom, processed_geom
                )
            
            metrics[theme_id] = theme_metrics
        
        return metrics
    
    def _prepare_properties_for_export(
        self,
        dict_processresults: Dict,
        include_formula: bool,
        include_attributes: bool
    ) -> Optional[Dict]:
        """
        Prepare properties dictionary for export.
        
        Args:
            dict_processresults: Processing results
            include_formula: Whether to include formula information
            include_attributes: Whether to include geometry attributes
            
        Returns:
            Properties dictionary or None
        """
        if not (include_formula or include_attributes):
            return None
        
        series_prop_dict = {}
        
        for theme_id, results in dict_processresults.items():
            theme_props = {}
            
            for distance, result in results.items():
                props = {}
                
                if include_formula and 'formula' in result:
                    props.update(result['formula'])
                
                if include_attributes and 'result' in result:
                    geom = result['result']
                    if hasattr(geom, 'area'):
                        props['area'] = geom.area
                    if hasattr(geom, 'length'):
                        props['perimeter'] = geom.length
                
                # Add any additional properties from the result
                for key, value in result.items():
                    if key not in ['result', 'result_diff', 'result_diff_plus', 'result_diff_min']:
                        if isinstance(value, (str, int, float, bool)):
                            props[key] = value
                
                theme_props[distance] = props
            
            series_prop_dict[theme_id] = theme_props
        
        return series_prop_dict
    
    def _calculate_single_metrics(self, original_geom, processed_geom) -> Dict:
        """
        Calculate metrics for a single geometry comparison.
        
        Args:
            original_geom: Original geometry
            processed_geom: Processed geometry
            
        Returns:
            Dictionary of calculated metrics
        """
        metrics = {}
        
        try:
            if hasattr(original_geom, 'area') and hasattr(processed_geom, 'area'):
                metrics['area_original'] = original_geom.area
                metrics['area_processed'] = processed_geom.area
                metrics['area_difference'] = abs(processed_geom.area - original_geom.area)
                
                if original_geom.area > 0:
                    metrics['area_change_percentage'] = (
                        (processed_geom.area - original_geom.area) / original_geom.area * 100
                    )
            
            if hasattr(original_geom, 'length') and hasattr(processed_geom, 'length'):
                metrics['perimeter_original'] = original_geom.length
                metrics['perimeter_processed'] = processed_geom.length
                metrics['perimeter_difference'] = abs(processed_geom.length - original_geom.length)
        
        except Exception as e:
            self.logger.feedback_warning(f"Error calculating metrics: {e}")
        
        return metrics
