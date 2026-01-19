# Comprehensive Plan to Fix Failing Tests

## Analysis of Test Failures

I've analyzed the failing tests and identified the root causes:

### 1. **Connection Errors in Multiple Tests**
**Tests affected:**
- `test_example_131635`
- `test_example_combined_borders_adp_gbg`
- `test_example_wanted_changes`
- `test_example_predictor`
- `test_load_data`
- `test_onroerenderfgoedloader_by_aanduidid`
- `test_onroerenderfgoedloader_by_erfgoedid`

**Root cause:** These tests are failing because they're trying to make HTTP requests that aren't being mocked properly. The error "Connection refused by Responses - the call doesn't match any registered mock" indicates that the test is trying to make a request that doesn't match any of the registered mock responses.

**Specific issues:**
- The tests are using `requests_mock` fixture but the URLs being called don't match the mock patterns
- Some tests are missing proper mock setup for the GRB API calls
- The `callback_grb_response` fixture is being used but may not cover all the required endpoints

### 2. **AttributeError in test_grbspecificdateparcelloader**
**Test affected:** `test_grbspecificdateparcelloader`

**Root cause:** In the `GRBSpecificDateParcelLoader` class, there's an issue where `date.year` is being called on what appears to be a string object instead of a date object. Looking at the code in `loader.py`, line 233:

```python
self.data_dict_source["source_url"] = GRB_FISCAL_PARCELS_URL + "/Adpf" + str(date.year)
```

The `date` parameter is expected to be a date object, but it's being treated as a string in some cases.

### 3. **Version Consistency Error**
**Test affected:** `test_version_consistency`

**Root cause:** The test is checking that the version in `pyproject.toml` matches the `__version__` in the package, but there's an assertion error indicating that not all requests have been executed. This suggests there might be some issue with the test setup or the version strings don't match.

## Comprehensive Plan to Fix All Test Failures

### Step 1: Fix the Mock Response Issues

**Files to modify:**
- `tests/test_examples.py`
- `tests/test_loader.py`
- `tests/test_oe.py`
- `tests/conftest.py`

**Actions:**
1. **Update the `callback_grb_response` fixture** to handle more URL patterns and ensure it covers all GRB API endpoints being called in the tests
2. **Add proper mock responses** for all the failing tests that make HTTP requests
3. **Ensure URL patterns match** exactly what the tests are calling
4. **Add missing mock responses** for the OnroerendErfgoedLoader tests

### Step 2: Fix the AttributeError in GRBSpecificDateParcelLoader

**File to modify:**
- `brdr/be/grb/loader.py`

**Actions:**
1. **Fix the date handling** in the `GRBSpecificDateParcelLoader.__init__` method
2. **Ensure the date parameter is properly converted** to a date object before accessing `.year`
3. **Add proper validation** to ensure the date parameter is in the correct format

### Step 3: Fix the Version Consistency Test

**File to modify:**
- `tests/test_version.py`

**Actions:**
1. **Debug the version consistency issue** by checking what versions are being compared
2. **Ensure the test is properly set up** and not making any unexpected requests
3. **Fix any version mismatch** between `pyproject.toml` and the package `__version__`

### Step 4: Run Tests to Verify Fixes

**Actions:**
1. **Run the specific failing tests** to verify each fix works
2. **Run the full test suite** to ensure no regressions are introduced
3. **Iterate on any remaining issues** that surface during testing

## Detailed Implementation Plan

### Fix 1: Update Mock Responses and URL Patterns

```python
# In tests/conftest.py, update the callback_grb_response fixture
@pytest.fixture
def callback_grb_response(requests_mock):
    response = copy(grb_responses.grb_response)

    def callback(request):
        json_response = json.dumps(response)
        return 200, {}, json_response

    # Add more URL patterns to cover all GRB API endpoints
    grb_urls = [
        "https://geo.api.vlaanderen.be",
        "https://geo.api.vlaanderen.be/GRB",
        "https://geo.api.vlaanderen.be/GRB/ogc",
        "https://geo.api.vlaanderen.be/GRB/wfs",
        "https://www.mercator.vlaanderen.be"
    ]

    requests_mock.add_callback(
        method=requests_mock.GET,
        url=multi_url_pattern(*grb_urls),
        callback=callback,
    )
    return response
```

### Fix 2: Update Failing Tests with Proper Mocks

For each failing test, add the appropriate mock responses:

```python
# Example for test_example_131635
def test_example_131635(self, requests_mock):
    # Add mock for Mercator WFS
    requests_mock.add(
        requests_mock.GET,
        "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature&TYPENAMES=ps%3Aps_aandobj&SRSNAME=http%3A%2F%2Fwww.opengis.net%2Fdef%2Fcrs%2FEPSG%2F0%2F31370&outputFormat=application%2Fjson&limit=10000&CQL_FILTER=aanduid_id+IN+%28131635%29",
        json=mercator_responses.response1,
        status=200,
    )

    # Add mock for GRB API calls
    requests_mock.add(
        requests_mock.GET,
        "https://geo.api.vlaanderen.be/GRB/ogc/features/v1/collections/ADP/items",
        json=grb_responses.grb_response,
        status=200,
    )

    # Rest of the test...
```

### Fix 3: Fix the Date Handling in GRBSpecificDateParcelLoader

```python
# In brdr/be/grb/loader.py
def __init__(self, date: str, aligner: Any, partition: int = 1000):
    # ... existing code ...

    # Fix the date handling
    if isinstance(date, str):
        try:
            date_obj = datetime.strptime(date, DATE_FORMAT).date()
        except ValueError:
            raise ValueError(
                f"No valid date, please provide a date in the format: {DATE_FORMAT}"
            )
    elif isinstance(date, date):
        date_obj = date
    else:
        raise ValueError(
            f"Date must be a string in format {DATE_FORMAT} or a date object"
        )

    # Use date_obj instead of date for year access
    self.data_dict_source["source_url"] = GRB_FISCAL_PARCELS_URL + "/Adpf" + str(date_obj.year)

    # ... rest of the code ...
```

### Fix 4: Debug and Fix Version Consistency

```python
# In tests/test_version.py
def test_version_consistency(pyproject_version):
    """Test to ensure version in pyproject.toml matches __version__."""
    print(f"Pyproject version: {pyproject_version}")
    print(f"Package version: {__version__}")
    assert pyproject_version == __version__
```

## Expected Outcomes

1. **All connection errors will be resolved** by properly mocking the HTTP requests
2. **The AttributeError will be fixed** by proper date handling
3. **Version consistency will be ensured** by debugging the version mismatch
4. **All tests will pass** after implementing these fixes

## Risk Assessment

- **Low risk** for the mock response fixes - these are isolated to test files
- **Medium risk** for the date handling fix - need to ensure backward compatibility
- **Low risk** for the version consistency fix - this is a simple comparison

## Timeline Estimate

- **Mock response fixes**: 1-2 hours
- **Date handling fix**: 30 minutes
- **Version consistency fix**: 30 minutes
- **Testing and verification**: 1-2 hours

**Total estimated time**: 3-5 hours

## Implementation Order

1. **Fix the mock response issues first** - this will resolve the majority of failing tests
2. **Fix the date handling issue** - this is a straightforward fix in the loader
3. **Fix the version consistency** - this should be quick to debug and fix
4. **Run comprehensive tests** to ensure everything works together

## Verification Strategy

After implementing each fix:
1. Run the specific test that was failing to verify the fix
2. Run related tests to ensure no regressions
3. After all fixes, run the full test suite
4. Document any edge cases or additional fixes needed