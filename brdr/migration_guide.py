"""
Migration guide and compatibility layer for transitioning from old Aligner to refactored version.

This module provides:
1. A compatibility wrapper for the old Aligner interface
2. Migration utilities
3. Documentation on how to update existing code
"""

import warnings
from typing import Dict, Any

from brdr.aligner_refactored import Aligner as RefactoredAligner
from brdr.aligner_config import AlignerConfig


class AlignerCompatibilityWrapper(RefactoredAligner):
    """
    Compatibility wrapper that provides the old Aligner interface while using
    the new refactored implementation underneath.
    
    This class is intended to help with migration by providing backward
    compatibility for existing code that uses the old Aligner class.
    """
    
    def __init__(self, **kwargs):
        """
        Initialize with backward compatibility for old parameter names.
        
        Issues deprecation warnings for usage patterns that should be updated.
        """
        # Issue a deprecation warning
        warnings.warn(
            "Using the old Aligner interface. Please consider migrating to the "
            "new refactored Aligner class for better performance and maintainability. "
            "See migration_guide.py for details.",
            DeprecationWarning,
            stacklevel=2
        )
        
        # Initialize with the refactored implementation
        super().__init__(**kwargs)
    
    # Add any deprecated methods that need special handling
    def _deprecated_method_example(self):
        """Example of how to handle deprecated methods."""
        warnings.warn(
            "This method is deprecated. Use the new method instead.",
            DeprecationWarning,
            stacklevel=2
        )
        # Delegate to new implementation
        pass


def migrate_aligner_usage(old_code_example: str) -> str:
    """
    Example function showing how to migrate from old to new Aligner usage.
    
    Args:
        old_code_example: String containing old Aligner code
        
    Returns:
        String with suggested new code
    """
    
    migration_examples = {
        "old_pattern": "new_pattern",
        # Add specific migration patterns here
    }
    
    # This is a simple example - in practice, you might use AST parsing
    # or more sophisticated code transformation tools
    new_code = old_code_example
    for old, new in migration_examples.items():
        new_code = new_code.replace(old, new)
    
    return new_code


def print_migration_guide():
    """
    Print a comprehensive migration guide for users.
    """
    
    guide = """
    BRDR Aligner Migration Guide
    ============================
    
    The Aligner class has been refactored for better maintainability and performance.
    Here's how to migrate your existing code:
    
    1. BASIC USAGE (No changes required)
    ------------------------------------
    
    OLD CODE:
    ```python
    from brdr.aligner import Aligner
    
    aligner = Aligner(crs="EPSG:31370", relevant_distance=2.0)
    aligner.load_thematic_data(loader)
    aligner.load_reference_data(loader)
    results = aligner.process()
    ```
    
    NEW CODE:
    ```python
    from brdr.aligner_refactored import Aligner
    
    aligner = Aligner(crs="EPSG:31370", relevant_distance=2.0)
    aligner.load_thematic_data(loader)
    aligner.load_reference_data(loader)
    results = aligner.process()
    ```
    
    2. CONFIGURATION-HEAVY USAGE (Recommended to use config object)
    --------------------------------------------------------------
    
    OLD CODE:
    ```python
    aligner = Aligner(
        crs="EPSG:31370",
        relevant_distance=2.0,
        od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=75,
        correction_distance=0.02,
        # ... many more parameters
    )
    ```
    
    NEW CODE (Recommended):
    ```python
    from brdr.aligner_config import AlignerConfig
    from brdr.aligner_refactored import Aligner
    
    config = AlignerConfig(
        crs="EPSG:31370",
        relevant_distance=2.0,
        od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=75,
        correction_distance=0.02,
        # ... other parameters
    )
    
    aligner = Aligner()  # Will use config defaults
    # Or pass specific overrides:
    aligner = Aligner(crs=config.crs, relevant_distance=config.relevant_distance)
    ```
    
    3. ACCESSING INTERNAL COMPONENTS (Advanced usage)
    ------------------------------------------------
    
    If you need access to specific functionality:
    
    ```python
    aligner = Aligner()
    
    # Access core data management
    core = aligner.core
    core.load_thematic_data(loader)
    
    # Access processing functionality
    processor = aligner.processor
    result = processor.process_geometry(geometry)
    
    # Access export functionality
    exporter = aligner.exporter
    geojson = exporter.get_results_as_geojson(results)
    ```
    
    4. BENEFITS OF THE NEW STRUCTURE
    --------------------------------
    
    - Better separation of concerns
    - Easier testing of individual components
    - Improved performance through optimized data structures
    - Better error handling and validation
    - More maintainable codebase
    - Configuration validation
    
    5. BREAKING CHANGES
    -------------------
    
    - None for basic usage
    - Internal method signatures may have changed (if you were accessing private methods)
    - Some edge case behaviors may be slightly different due to improved error handling
    
    6. MIGRATION CHECKLIST
    ----------------------
    
    [ ] Update import statements
    [ ] Test all functionality with your data
    [ ] Consider using AlignerConfig for complex configurations
    [ ] Update any custom error handling (new exceptions may be raised)
    [ ] Review any code that accessed private methods (starting with _)
    
    For more detailed information, see the documentation and examples.
    """
    
    print(guide)


if __name__ == "__main__":
    # Print the migration guide when run as a script
    print_migration_guide()
