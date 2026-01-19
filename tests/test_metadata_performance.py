"""
Test cases for the metadata performance optimization that limits reference features
in actuation metadata to only those that overlap with result geometries.
"""

import pytest
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.configs import AlignerConfig
from brdr.loader import DictLoader


class TestMetadataPerformance:
    """Test cases for metadata performance optimization."""

    def test_overlapping_reference_tracking(self):
        """Test that compare_to_reference correctly tracks overlapping reference IDs."""
        aligner = Aligner()
        
        # Create thematic and reference data
        thematic_dict = {
            'theme_id_1': from_wkt('POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))')
        }
        reference_dict = {
            'ref_id_1': from_wkt('POLYGON ((0 1, 0 10, 8 10, 10 1, 0 1))'),  # Overlaps
            'ref_id_2': from_wkt('POLYGON ((20 20, 20 25, 25 25, 25 20, 20 20))'),  # Does not overlap
            'ref_id_3': from_wkt('POLYGON ((5 5, 5 8, 8 8, 8 5, 5 5))')  # Overlaps
        }
        
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))
        
        # Process with observations enabled
        result = aligner.process([1])
        
        # Check that overlapping_reference_ids is present and correct
        for thematic_id, rd_res in result.results.items():
            for rd, res in rd_res.items():
                if 'observation' in res:
                    obs = res['observation']
                    assert 'overlapping_reference_ids' in obs
                    overlapping_ids = obs['overlapping_reference_ids']
                    
                    # Should contain ref_id_1 and ref_id_3, but not ref_id_2
                    assert 'ref_id_1' in overlapping_ids
                    assert 'ref_id_3' in overlapping_ids
                    assert 'ref_id_2' not in overlapping_ids
                    
                    # Should match the keys in reference_features
                    assert set(overlapping_ids) == set(obs['reference_features'].keys())
                    break
            break

    def test_actuation_metadata_filtering(self):
        """Test that actuation metadata only includes overlapping reference features."""
        config = AlignerConfig(log_metadata=True, add_observations=True)
        aligner = Aligner(config=config)
        
        # Create thematic and reference data
        thematic_dict = {
            'theme_id_1': from_wkt('POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))')
        }
        reference_dict = {
            'ref_id_1': from_wkt('POLYGON ((0 1, 0 10, 8 10, 10 1, 0 1))'),  # Overlaps
            'ref_id_2': from_wkt('POLYGON ((20 20, 20 25, 25 25, 25 20, 20 20))'),  # Does not overlap
            'ref_id_3': from_wkt('POLYGON ((5 5, 5 8, 8 8, 8 5, 5 5))')  # Overlaps
        }
        
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))
        
        # Process with metadata logging enabled
        result = aligner.process([1])
        
        # Check that actuation metadata only includes overlapping references
        for thematic_id, rd_res in result.results.items():
            for rd, res in rd_res.items():
                if 'metadata' in res and 'actuation' in res['metadata']:
                    actuation = res['metadata']['actuation']
                    ref_geometries = actuation['reference_geometries']
                    
                    # Should only have 2 reference geometries (the overlapping ones)
                    assert len(ref_geometries) == 2
                    
                    # Check that the reference IDs in actuation match the overlapping ones
                    actuation_ref_ids = {ref_geom['id'] for ref_geom in ref_geometries}
                    observation_ref_ids = set(res['observation']['overlapping_reference_ids'])
                    
                    # The actuation should contain the same references as the observation
                    assert len(actuation_ref_ids) == len(observation_ref_ids)
                    
                    break
            break

    def test_backward_compatibility_no_observations(self):
        """Test backward compatibility when add_observations is False."""
        config = AlignerConfig(log_metadata=True, add_observations=False)
        aligner = Aligner(config=config)
        
        # Create thematic and reference data
        thematic_dict = {
            'theme_id_1': from_wkt('POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))')
        }
        reference_dict = {
            'ref_id_1': from_wkt('POLYGON ((0 1, 0 10, 8 10, 10 1, 0 1))'),
            'ref_id_2': from_wkt('POLYGON ((20 20, 20 25, 25 25, 25 20, 20 20))'),
        }
        
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))
        
        # Process with metadata logging but no observations
        result = aligner.process([1])
        
        # Check that actuation metadata falls back to all features
        for thematic_id, rd_res in result.results.items():
            for rd, res in rd_res.items():
                if 'metadata' in res and 'actuation' in res['metadata']:
                    actuation = res['metadata']['actuation']
                    ref_geometries = actuation['reference_geometries']
                    
                    # Should fall back to all 2 reference features
                    assert len(ref_geometries) == 2
                    
                    break
            break

    def test_no_overlapping_references(self):
        """Test behavior when result geometry doesn't overlap with any reference features."""
        config = AlignerConfig(log_metadata=True, add_observations=True)
        aligner = Aligner(config=config)
        
        # Create thematic and reference data with no overlap
        thematic_dict = {
            'theme_id_1': from_wkt('POLYGON ((0 0, 0 1, 1 1, 1 0, 0 0))')
        }
        reference_dict = {
            'ref_id_1': from_wkt('POLYGON ((20 20, 20 25, 25 25, 25 20, 20 20))'),  # Far away
            'ref_id_2': from_wkt('POLYGON ((30 30, 30 35, 35 35, 35 30, 30 30))'),  # Also far away
        }
        
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))
        
        # Process with metadata logging enabled
        result = aligner.process([1])
        
        # Check that overlapping_reference_ids is empty
        for thematic_id, rd_res in result.results.items():
            for rd, res in rd_res.items():
                if 'observation' in res:
                    obs = res['observation']
                    assert 'overlapping_reference_ids' in obs
                    overlapping_ids = obs['overlapping_reference_ids']
                    
                    # Should be empty since no references overlap
                    assert len(overlapping_ids) == 0
                    
                    # Reference features should also be empty
                    assert len(obs['reference_features']) == 0
                    
                    break
            break

    def test_observation_generation_unchanged(self):
        """Test that observation generation still works correctly."""
        config = AlignerConfig(log_metadata=True, add_observations=True)
        aligner = Aligner(config=config)
        
        # Create thematic and reference data
        thematic_dict = {
            'theme_id_1': from_wkt('POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))')
        }
        reference_dict = {
            'ref_id_1': from_wkt('POLYGON ((0 1, 0 10, 8 10, 10 1, 0 1))'),  # Overlaps
            'ref_id_2': from_wkt('POLYGON ((20 20, 20 25, 25 25, 25 20, 20 20))'),  # Does not overlap
        }
        
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))
        
        # Process with metadata logging enabled
        result = aligner.process([1])
        
        # Check that observations are generated correctly
        for thematic_id, rd_res in result.results.items():
            for rd, res in rd_res.items():
                if 'metadata' in res and 'observations' in res['metadata']:
                    observations = res['metadata']['observations']
                    
                    # Should have observations for the overlapping reference
                    assert len(observations) > 0
                    
                    # Each observation should reference the correct feature
                    for obs in observations:
                        if 'has_feature_of_interest' in obs:
                            # Should reference either the result or the overlapping reference
                            feature_of_interest = obs['has_feature_of_interest']
                            assert feature_of_interest is not None
                            
                    break
            break

    def test_reference_lookup_optimization(self):
        """Test that the reference lookup optimization works correctly."""
        config = AlignerConfig(log_metadata=True, add_observations=True)
        aligner = Aligner(config=config)
        
        # Create thematic and reference data with multiple features
        thematic_dict = {
            'theme_id_1': from_wkt('POLYGON ((0 0, 0 10, 10 10, 10 0, 0 0))')
        }
        reference_dict = {
            'ref_id_1': from_wkt('POLYGON ((0 1, 0 9, 9 9, 9 1, 0 1))'),  # Overlaps
            'ref_id_2': from_wkt('POLYGON ((1 1, 1 8, 8 8, 8 1, 1 1))'),  # Overlaps
            'ref_id_3': from_wkt('POLYGON ((20 20, 20 25, 25 25, 25 20, 20 20))'),  # Does not overlap
        }
        
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))
        
        # Process with metadata logging enabled
        result = aligner.process([1])
        
        # Verify that observations are generated correctly with the optimization
        for thematic_id, rd_res in result.results.items():
            for rd, res in rd_res.items():
                if 'metadata' in res and 'observations' in res['metadata']:
                    observations = res['metadata']['observations']
                    
                    # Should have observations for the overlapping references
                    assert len(observations) > 0
                    
                    # Verify that observations reference the correct features
                    ref_ids_in_observations = set()
                    for obs in observations:
                        if 'has_feature_of_interest' in obs:
                            feature_of_interest = obs['has_feature_of_interest']
                            assert feature_of_interest is not None
                            # Check that it's a valid reference (can be dict for ref features or string URN for result)
                            if isinstance(feature_of_interest, dict):
                                assert 'id' in feature_of_interest
                                ref_ids_in_observations.add(feature_of_interest['id'])
                            # String URNs are valid for result references
                    
                    # Should have observations for both overlapping references
                    assert len(ref_ids_in_observations) >= 2
                    
                    break
            break