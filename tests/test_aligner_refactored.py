"""
Tests for the refactored Aligner implementation.
"""
import unittest
from unittest.mock import Mock, patch

from shapely.geometry import Polygon
from shapely import from_wkt

from brdr.aligner_config import AlignerConfig
from brdr.aligner_core import AlignerCore
from brdr.aligner_processor import AlignerProcessor
from brdr.aligner_exporter import AlignerExporter
from brdr.aligner_refactored import Aligner
from brdr.loader import DictLoader
from brdr.enums import OpenDomainStrategy, AlignerInputType


class TestAlignerConfig(unittest.TestCase):
    """Test the AlignerConfig class."""
    
    def test_default_config(self):
        """Test default configuration values."""
        config = AlignerConfig()
        
        self.assertEqual(config.relevant_distance, 1.0)
        self.assertEqual(config.crs, "EPSG:31370")
        self.assertEqual(config.od_strategy, OpenDomainStrategy.SNAP_ALL_SIDE)
        self.assertEqual(config.threshold_overlap_percentage, 50)
        self.assertTrue(config.multi_as_single_modus)
        self.assertFalse(config.preserve_topology)
    
    def test_config_validation(self):
        """Test configuration validation."""
        # Test invalid relevant_distance
        with self.assertRaises(ValueError):
            AlignerConfig(relevant_distance=-1)
        
        # Test invalid threshold_overlap_percentage
        with self.assertRaises(ValueError):
            AlignerConfig(threshold_overlap_percentage=150)
        
        # Test invalid threshold_circle_ratio
        with self.assertRaises(ValueError):
            AlignerConfig(threshold_circle_ratio=1.5)
    
    def test_custom_config(self):
        """Test custom configuration values."""
        config = AlignerConfig(
            relevant_distance=2.5,
            crs="EPSG:4326",
            threshold_overlap_percentage=75,
            correction_distance=0.02
        )
        
        self.assertEqual(config.relevant_distance, 2.5)
        self.assertEqual(config.crs, "EPSG:4326")
        self.assertEqual(config.threshold_overlap_percentage, 75)
        self.assertEqual(config.correction_distance, 0.02)


class TestAlignerCore(unittest.TestCase):
    """Test the AlignerCore class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.config = AlignerConfig()
        self.core = AlignerCore(self.config)
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
    
    def test_initialization(self):
        """Test core initialization."""
        self.assertEqual(self.core.config, self.config)
        self.assertEqual(len(self.core.dict_thematic), 0)
        self.assertEqual(len(self.core.dict_reference), 0)
        self.assertIsNone(self.core.thematic_union)
        self.assertIsNone(self.core.reference_union)
    
    def test_load_thematic_data(self):
        """Test loading thematic data."""
        thematic_dict = {"theme_1": self.sample_geom}
        loader = DictLoader(thematic_dict)
        
        self.core.load_thematic_data(loader)
        
        self.assertEqual(len(self.core.dict_thematic), 1)
        self.assertIn("theme_1", self.core.dict_thematic)
        self.assertEqual(self.core.dict_thematic["theme_1"], self.sample_geom)
    
    def test_load_reference_data(self):
        """Test loading reference data."""
        reference_dict = {"ref_1": self.sample_geom}
        loader = DictLoader(reference_dict)
        
        self.core.load_reference_data(loader)
        
        self.assertEqual(len(self.core.dict_reference), 1)
        self.assertIn("ref_1", self.core.dict_reference)
        self.assertEqual(self.core.dict_reference["ref_1"], self.sample_geom)
        self.assertIsNotNone(self.core.reference_tree)
        self.assertIsNotNone(self.core.reference_items)
    
    def test_load_empty_data(self):
        """Test loading empty data raises ValueError."""
        empty_loader = DictLoader({})
        
        with self.assertRaises(ValueError):
            self.core.load_thematic_data(empty_loader)
        
        with self.assertRaises(ValueError):
            self.core.load_reference_data(empty_loader)
    
    def test_load_none_loader(self):
        """Test loading with None loader raises ValueError."""
        with self.assertRaises(ValueError):
            self.core.load_thematic_data(None)
        
        with self.assertRaises(ValueError):
            self.core.load_reference_data(None)
    
    def test_get_thematic_union(self):
        """Test getting thematic union."""
        thematic_dict = {"theme_1": self.sample_geom}
        loader = DictLoader(thematic_dict)
        self.core.load_thematic_data(loader)
        
        union = self.core.get_thematic_union()
        
        self.assertIsNotNone(union)
        self.assertEqual(union, self.sample_geom)
    
    def test_get_thematic_union_no_data(self):
        """Test getting thematic union with no data raises ValueError."""
        with self.assertRaises(ValueError):
            self.core.get_thematic_union()
    
    def test_validate_data_loaded(self):
        """Test data validation."""
        # No data loaded
        with self.assertRaises(ValueError):
            self.core.validate_data_loaded()
        
        # Only thematic data loaded
        thematic_dict = {"theme_1": self.sample_geom}
        loader = DictLoader(thematic_dict)
        self.core.load_thematic_data(loader)
        
        with self.assertRaises(ValueError):
            self.core.validate_data_loaded()
        
        # Both data types loaded
        reference_dict = {"ref_1": self.sample_geom}
        loader = DictLoader(reference_dict)
        self.core.load_reference_data(loader)
        
        # Should not raise
        self.core.validate_data_loaded()


class TestAlignerRefactored(unittest.TestCase):
    """Test the refactored Aligner class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
    
    def test_initialization(self):
        """Test aligner initialization."""
        self.assertIsInstance(self.aligner.config, AlignerConfig)
        self.assertIsInstance(self.aligner.core, AlignerCore)
        self.assertIsInstance(self.aligner.processor, AlignerProcessor)
        self.assertIsInstance(self.aligner.exporter, AlignerExporter)
    
    def test_initialization_with_parameters(self):
        """Test aligner initialization with custom parameters."""
        aligner = Aligner(
            crs="EPSG:4326",
            relevant_distance=2.5,
            threshold_overlap_percentage=75
        )
        
        self.assertEqual(aligner.config.crs, "EPSG:4326")
        self.assertEqual(aligner.config.relevant_distance, 2.5)
        self.assertEqual(aligner.config.threshold_overlap_percentage, 75)
    
    def test_load_data_methods(self):
        """Test data loading methods."""
        thematic_dict = {"theme_1": self.sample_geom}
        reference_dict = {"ref_1": self.sample_geom}
        
        thematic_loader = DictLoader(thematic_dict)
        reference_loader = DictLoader(reference_dict)
        
        self.aligner.load_thematic_data(thematic_loader)
        self.aligner.load_reference_data(reference_loader)
        
        self.assertEqual(len(self.aligner.dict_thematic), 1)
        self.assertEqual(len(self.aligner.dict_reference), 1)
    
    def test_properties(self):
        """Test property access."""
        # Test CRS property
        self.assertEqual(self.aligner.CRS, "EPSG:31370")
        
        # Test data properties
        self.assertEqual(len(self.aligner.dict_thematic), 0)
        self.assertEqual(len(self.aligner.dict_reference), 0)
        
        # Test configuration properties
        self.assertEqual(self.aligner.relevant_distances, self.aligner.config.relevant_distances)
        self.assertEqual(self.aligner.od_strategy, self.aligner.config.od_strategy)
        self.assertEqual(self.aligner.threshold_overlap_percentage, self.aligner.config.threshold_overlap_percentage)
    
    def test_property_setters(self):
        """Test property setters."""
        # Test setting relevant_distances
        new_distances = [0.5, 1.0, 1.5]
        self.aligner.relevant_distances = new_distances
        self.assertEqual(self.aligner.relevant_distances, new_distances)
        
        # Test setting od_strategy
        new_strategy = OpenDomainStrategy.EXCLUDE
        self.aligner.od_strategy = new_strategy
        self.assertEqual(self.aligner.od_strategy, new_strategy)
        
        # Test setting threshold_overlap_percentage
        new_threshold = 75
        self.aligner.threshold_overlap_percentage = new_threshold
        self.assertEqual(self.aligner.threshold_overlap_percentage, new_threshold)
    
    def test_get_input_as_geojson(self):
        """Test getting input as GeoJSON."""
        # Load some test data
        thematic_dict = {"theme_1": self.sample_geom}
        reference_dict = {"ref_1": self.sample_geom}
        
        self.aligner.load_thematic_data(DictLoader(thematic_dict))
        self.aligner.load_reference_data(DictLoader(reference_dict))
        
        # Test getting reference data
        ref_geojson = self.aligner.get_input_as_geojson(AlignerInputType.REFERENCE)
        self.assertIsNotNone(ref_geojson)
        
        # Test getting thematic data
        theme_geojson = self.aligner.get_input_as_geojson(AlignerInputType.THEMATIC)
        self.assertIsNotNone(theme_geojson)
    
    def test_repr(self):
        """Test string representation."""
        repr_str = repr(self.aligner)
        self.assertIn("Aligner", repr_str)
        self.assertIn("thematic_features=0", repr_str)
        self.assertIn("reference_features=0", repr_str)
        self.assertIn("crs=EPSG:31370", repr_str)


class TestBackwardCompatibility(unittest.TestCase):
    """Test backward compatibility with the original Aligner interface."""
    
    def test_basic_usage_compatibility(self):
        """Test that basic usage patterns still work."""
        # This test ensures that existing code will continue to work
        aligner = Aligner(
            crs="EPSG:31370",
            relevant_distance=2.0,
            od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE
        )
        
        sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        thematic_dict = {"theme_1": sample_geom}
        reference_dict = {"ref_1": sample_geom}
        
        # Load data
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))
        
        # Access properties
        self.assertEqual(aligner.CRS, "EPSG:31370")
        self.assertEqual(len(aligner.dict_thematic), 1)
        self.assertEqual(len(aligner.dict_reference), 1)
        
        # Get union
        union = aligner.get_thematic_union()
        self.assertIsNotNone(union)


if __name__ == "__main__":
    unittest.main()

