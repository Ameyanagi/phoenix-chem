# Batch Processing

PHOENIX provides efficient batch screening for analyzing multiple compounds.

## Basic Usage

```python
from phoenix import screen

smiles_list = [
    "CCO",                                    # Ethanol
    "CC(=O)C",                                # Acetone
    "c1ccccc1[N+](=O)[O-]",                   # Nitrobenzene
    "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",  # TNT
]

results = screen(smiles_list)
print(f"Screened {results.successful} compounds successfully")
print(f"Failed: {results.failed}")
```

## BatchResult

The `screen()` function returns a `BatchResult` object:

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `dataframe` | pd.DataFrame | Results DataFrame |
| `successful` | int | Number of successful screenings |
| `failed` | int | Number of failed screenings |

### Methods

| Method | Description |
|--------|-------------|
| `to_csv(path)` | Export to CSV file |
| `to_json()` | Export to JSON string |
| `summary()` | Get summary statistics |

## DataFrame Columns

The results DataFrame contains:

| Column | Type | Description |
|--------|------|-------------|
| `smiles` | str | Input SMILES |
| `canonical_smiles` | str | Canonicalized SMILES |
| `formula` | str | Molecular formula |
| `mw` | float | Molecular weight (g/mol) |
| `delta_hf_kJ_mol` | float | ΔHf° (kJ/mol) |
| `delta_hd_kJ_mol` | float | ΔHd (kJ/mol) |
| `delta_hd_cal_g` | float | ΔHd (cal/g) |
| `ob_percent` | float | Oxygen balance (%) |
| `hazard_class` | str | HIGH, MEDIUM, or LOW |
| `triggered_criteria` | list | CHETAH criteria triggered |
| `alerts` | list | Functional group alerts |
| `gas_volume_L_g` | float | Gas generation (L/g) |
| `error` | str | Error type if failed |
| `error_message` | str | Error details if failed |

## Working with Results

### Accessing the DataFrame

```python
results = screen(smiles_list)
df = results.dataframe

# View all columns
print(df.columns.tolist())

# Basic statistics
print(df.describe())
```

### Filtering Results

```python
df = results.dataframe

# High hazard compounds
high_hazard = df[df['hazard_class'] == 'HIGH']
print(f"High hazard: {len(high_hazard)}")

# Compounds with specific criteria
criterion_1 = df[df['triggered_criteria'].apply(lambda x: 1 in x if x else False)]

# Compounds with alerts
has_alerts = df[df['alerts'].apply(lambda x: len(x) > 0 if x else False)]
```

### Sorting Results

```python
df = results.dataframe

# Sort by ΔHd (most exothermic first)
sorted_df = df.sort_values('delta_hd_cal_g')

# Sort by hazard class
hazard_order = {'HIGH': 0, 'MEDIUM': 1, 'LOW': 2, None: 3}
df['hazard_rank'] = df['hazard_class'].map(hazard_order)
sorted_df = df.sort_values('hazard_rank')
```

## Exporting Results

### To CSV

```python
results = screen(smiles_list)
results.to_csv("screening_results.csv")
```

### To JSON

```python
json_str = results.to_json()
print(json_str)

# Write to file
with open("results.json", "w") as f:
    f.write(json_str)
```

### To Excel

```python
df = results.dataframe
df.to_excel("results.xlsx", index=False)
```

## Summary Statistics

```python
results = screen(smiles_list)
summary = results.summary()

print(f"Total compounds: {summary['total_compounds']}")
print(f"Successful: {summary['successful']}")
print(f"Failed: {summary['failed']}")
print(f"Hazard distribution: {summary['hazard_class_counts']}")
print(f"ΔHd range: {summary['delta_hd_cal_g_min']:.1f} to {summary['delta_hd_cal_g_max']:.1f} cal/g")
```

## Progress Reporting

Use a callback for progress updates:

```python
def progress_callback(current, total):
    percent = (current / total) * 100
    print(f"\rProgress: {current}/{total} ({percent:.1f}%)", end="")

results = screen(smiles_list, progress_callback=progress_callback)
print()  # Newline after progress
```

### With tqdm

```python
from tqdm import tqdm

pbar = tqdm(total=len(smiles_list))

def update_progress(current, total):
    pbar.update(1)

results = screen(smiles_list, progress_callback=update_progress)
pbar.close()
```

## Error Handling

Failed compounds are captured in the DataFrame:

```python
results = screen(smiles_list)
df = results.dataframe

# Find failed compounds
failed = df[df['error'].notna()]

for _, row in failed.iterrows():
    print(f"SMILES: {row['smiles']}")
    print(f"Error: {row['error']}")
    print(f"Message: {row['error_message']}")
```

### Common Errors

| Error | Cause |
|-------|-------|
| `InvalidSmilesError` | Invalid SMILES string |
| `UnsupportedElementError` | Contains unsupported element |
| `UnsupportedStructureError` | Charged/radical species |
| `MissingGroupError` | Benson GA data unavailable |

### Graceful Degradation

Batch processing continues even when individual compounds fail:

```python
smiles_list = [
    "CCO",           # Valid
    "invalid",       # Invalid SMILES
    "[Fe]",          # Unsupported element
    "c1ccccc1",      # Valid
]

results = screen(smiles_list)

print(f"Successful: {results.successful}")  # 2
print(f"Failed: {results.failed}")          # 2
```

## Large-Scale Screening

### Memory Considerations

For very large datasets, process in chunks:

```python
import pandas as pd

def screen_in_chunks(all_smiles, chunk_size=1000):
    all_dfs = []

    for i in range(0, len(all_smiles), chunk_size):
        chunk = all_smiles[i:i+chunk_size]
        results = screen(chunk)
        all_dfs.append(results.dataframe)
        print(f"Processed {min(i+chunk_size, len(all_smiles))}/{len(all_smiles)}")

    return pd.concat(all_dfs, ignore_index=True)

# Usage
all_smiles = [...]  # Large list
combined_df = screen_in_chunks(all_smiles)
```

### Parallel Processing

For CPU-bound screening:

```python
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

def screen_chunk(smiles_chunk):
    return screen(smiles_chunk).dataframe

def parallel_screen(all_smiles, n_workers=4, chunk_size=100):
    chunks = [all_smiles[i:i+chunk_size]
              for i in range(0, len(all_smiles), chunk_size)]

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        dfs = list(executor.map(screen_chunk, chunks))

    return pd.concat(dfs, ignore_index=True)

# Usage
combined_df = parallel_screen(all_smiles)
```

## Practical Examples

### Screening a Chemical Library

```python
from phoenix import screen

# Read SMILES from file
with open("compounds.smi") as f:
    smiles_list = [line.strip() for line in f if line.strip()]

# Screen all compounds
results = screen(smiles_list)

# Save results
results.to_csv("screening_results.csv")

# Print summary
summary = results.summary()
print(f"\nScreening Summary:")
print(f"  Total: {summary['total_compounds']}")
print(f"  Successful: {summary['successful']}")
print(f"  High hazard: {summary['hazard_class_counts'].get('HIGH', 0)}")
print(f"  Medium hazard: {summary['hazard_class_counts'].get('MEDIUM', 0)}")
print(f"  Low hazard: {summary['hazard_class_counts'].get('LOW', 0)}")
```

### Identifying High-Risk Compounds

```python
results = screen(smiles_list)
df = results.dataframe

# Find high-risk compounds
high_risk = df[
    (df['hazard_class'] == 'HIGH') |
    (df['delta_hd_cal_g'] < -200)
]

print("High-risk compounds:")
for _, row in high_risk.iterrows():
    print(f"  {row['formula']}: {row['hazard_class']} ({row['delta_hd_cal_g']:.0f} cal/g)")
```

### Hazard Statistics Report

```python
results = screen(smiles_list)
df = results.dataframe

# Generate report
print("Hazard Screening Report")
print("=" * 50)
print(f"Compounds screened: {len(df)}")
print(f"Successful: {results.successful}")
print(f"Failed: {results.failed}")

print("\nHazard Class Distribution:")
for hazard, count in df['hazard_class'].value_counts().items():
    pct = count / results.successful * 100
    print(f"  {hazard}: {count} ({pct:.1f}%)")

print("\nΔHd Statistics:")
successful_df = df[df['error'].isna()]
print(f"  Min: {successful_df['delta_hd_cal_g'].min():.1f} cal/g")
print(f"  Max: {successful_df['delta_hd_cal_g'].max():.1f} cal/g")
print(f"  Mean: {successful_df['delta_hd_cal_g'].mean():.1f} cal/g")

print("\nMost Common Alerts:")
all_alerts = []
for alerts in successful_df['alerts'].dropna():
    all_alerts.extend(alerts)
from collections import Counter
for alert, count in Counter(all_alerts).most_common(5):
    print(f"  {alert}: {count}")
```

---

## Next Steps

- [Error Handling](error-handling.md): Handling exceptions
- [Hazard Evaluation](hazard-evaluation.md): Individual compound analysis
- [API Reference](../api/batch.md): BatchResult API
