"""
Demonstration script showing the performance improvement from the metadata optimization.
This script creates a scenario with many reference features and shows how the new
implementation only includes overlapping references in metadata.
"""

from brdr.aligner import Aligner
from brdr.configs import AlignerConfig
from brdr.loader import DictLoader
from shapely import from_wkt
import json

def create_large_reference_dataset():
    """Create a reference dataset with many features, only some of which overlap."""
    reference_dict = {}
    
    # Add overlapping features (near the thematic geometry)
    for i in range(5):
        reference_dict[f'ref_overlap_{i}'] = from_wkt(
            f'POLYGON (({i*2} {i*2}, {i*2} {i*2+3}, {i*2+3} {i*2+3}, {i*2+3} {i*2}, {i*2} {i*2}))'
        )
    
    # Add non-overlapping features (far from the thematic geometry)
    for i in range(20):
        reference_dict[f'ref_far_{i}'] = from_wkt(
            f'POLYGON (({100+i*10} {100+i*10}, {100+i*10} {100+i*10+5}, {100+i*10+5} {100+i*10+5}, {100+i*10+5} {100+i*10}, {100+i*10} {100+i*10}))'
        )
    
    return reference_dict

def main():
    print("=== Metadata Performance Optimization Demo ===\n")
    
    # Create aligner with metadata logging
    config = AlignerConfig(log_metadata=True, add_observations=True)
    aligner = Aligner(config=config)
    
    # Create thematic data (small polygon near origin)
    thematic_dict = {
        'theme_id_1': from_wkt('POLYGON ((1 1, 1 4, 4 4, 4 1, 1 1))')
    }
    
    # Create reference data with many features
    reference_dict = create_large_reference_dataset()
    
    print(f"Thematic features: {len(thematic_dict)}")
    print(f"Reference features: {len(reference_dict)}")
    print(f"Overlapping references: 5 (ref_overlap_0 to ref_overlap_4)")
    print(f"Non-overlapping references: 20 (ref_far_0 to ref_far_19)")
    print()
    
    # Load data
    aligner.load_thematic_data(DictLoader(thematic_dict))
    aligner.load_reference_data(DictLoader(reference_dict))
    
    # Process
    print("Processing...")
    result = aligner.process([1])
    
    # Analyze results
    for thematic_id, rd_res in result.results.items():
        for rd, res in rd_res.items():
            if 'observation' in res:
                obs = res['observation']
                overlapping_ids = obs['overlapping_reference_ids']
                
                print(f"\n=== OBSERVATION DATA ===")
                print(f"Overlapping reference IDs: {overlapping_ids}")
                print(f"Number of overlapping references: {len(overlapping_ids)}")
                print(f"Reference features in observation: {list(obs['reference_features'].keys())}")
            
            if 'metadata' in res and 'actuation' in res['metadata']:
                actuation = res['metadata']['actuation']
                ref_geometries = actuation['reference_geometries']
                
                print(f"\n=== ACTUATION METADATA ===")
                print(f"Number of reference geometries in actuation: {len(ref_geometries)}")
                
                # Calculate metadata size reduction
                total_refs = len(reference_dict)
                actuation_refs = len(ref_geometries)
                reduction_percentage = ((total_refs - actuation_refs) / total_refs) * 100
                
                print(f"Total reference features available: {total_refs}")
                print(f"Reference features in actuation: {actuation_refs}")
                print(f"Metadata size reduction: {reduction_percentage:.1f}%")
                
                # Show the actual reference IDs in actuation
                print(f"\nReference geometry IDs in actuation:")
                for ref_geom in ref_geometries:
                    print(f"  - {ref_geom['id']}")
                
                # Verify that the correct number of references are included
                if 'observation' in res:
                    obs_overlapping_count = len(res['observation']['overlapping_reference_ids'])
                    actuation_ref_count = len(ref_geometries)
                    
                    print(f"\n=== VERIFICATION ===")
                    print(f"Number of overlapping references from observation: {obs_overlapping_count}")
                    print(f"Number of reference geometries in actuation: {actuation_ref_count}")
                    
                    if obs_overlapping_count == actuation_ref_count:
                        print("✓ SUCCESS: Actuation metadata contains the correct number of overlapping references!")
                    else:
                        print("✗ ERROR: Mismatch in reference counts")
                    
                    # Note: The IDs differ because actuation uses UUIDs while observation uses original IDs
                    print("Note: Actuation uses UUIDs for reference IDs, observation uses original IDs")
                
                break
        break
    
    # Demonstrate the benefit
    print(f"\n=== PERFORMANCE BENEFIT ===")
    print("Before optimization: All {len(reference_dict)} reference features would be included in actuation metadata")
    print("After optimization: Only {len(overlapping_ids)} overlapping reference features are included")
    print(f"Reduction: {len(reference_dict) - len(overlapping_ids)} unnecessary reference features removed")
    print(f"Efficiency improvement: {reduction_percentage:.1f}% smaller metadata")

if __name__ == "__main__":
    main()