# Fix FragmentIntensityPrediction Test Failures

## Problem
4 tests in `TestFragmentIntensityModelBase.cs` fail with `Assert.That(isValid, Is.True)` returning `False`:
- `TestCollisionEnergyValidation` (line 86)
- `TestInstrumentTypeValidation` (line 120)
- `TestFragmentationTypeValidation` (line 155)
- `TestPrecursorChargeValidation` (line 188)

## Root Cause
In `TestFragmentIntensityModel` constructor (lines 49-51), untested constraints default to empty sets:
```csharp
AllowedCollisionEnergies = allowedEnergies ?? new HashSet<int>();
AllowedInstrumentTypes = allowedInstruments ?? new HashSet<string>();
AllowedFragmentationTypes = allowedFragmentations ?? new HashSet<string>();
```

Per `FragmentIntensityModel.cs` documentation:
- `null` = skip validation (not applicable)
- `empty set` = parameter IS required but any value accepted
- `populated set` = only listed values accepted

When testing one constraint, others default to empty sets (meaning "required"), causing validation to fail when test inputs have `null` for those fields.

## Fix
Change `Test/KoinaTests/FragmentIntensityPrediction/TestFragmentIntensityModelBase.cs` lines 49-51:
```csharp
// Before:
AllowedCollisionEnergies = allowedEnergies ?? new HashSet<int>();
AllowedInstrumentTypes = allowedInstruments ?? new HashSet<string>();
AllowedFragmentationTypes = allowedFragmentations ?? new HashSet<string>();

// After:
AllowedCollisionEnergies = allowedEnergies;
AllowedInstrumentTypes = allowedInstruments;
AllowedFragmentationTypes = allowedFragmentations;
```

This ensures untested constraints are skipped (`null`), while tests that explicitly pass empty sets (like `TestValidation_EmptyConstraintSets_AllowsAnyValue`) still work correctly.

## Verification
Run: `dotnet test ./Test/Test.csproj --filter "FullyQualifiedName~FragmentIntensityPrediction" --verbosity normal`
