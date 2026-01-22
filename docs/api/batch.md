# Batch API

Batch screening for multiple compounds.

---

## screen

::: phoenix.batch.screening.screen
    options:
      show_root_heading: true
      show_source: false

```python
from phoenix import screen

def screen(
    smiles_list: list[str],
    progress_callback: Callable[[int, int], None] | None = None,
) -> BatchResult
```

Screen multiple compounds for reactive hazards.

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `smiles_list` | list[str] | required | SMILES strings to process |
| `progress_callback` | Callable | None | Progress reporter |

### Returns

`BatchResult` - Results with DataFrame and statistics.

### Example

```python
from phoenix import screen

smiles_list = ["CCO", "c1ccccc1[N+](=O)[O-]", "invalid"]
results = screen(smiles_list)

print(f"Successful: {results.successful}")
print(f"Failed: {results.failed}")
```

---

## BatchResult

::: phoenix.batch.screening.BatchResult
    options:
      show_root_heading: true
      show_source: false

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `dataframe` | pd.DataFrame | Results DataFrame |
| `successful` | int | Successful screenings |
| `failed` | int | Failed screenings |

### Methods

#### to_csv

```python
def to_csv(self, path: str, **kwargs) -> None
```

Export results to CSV file.

#### to_json

```python
def to_json(self, **kwargs) -> str
```

Export results to JSON string.

#### summary

```python
def summary(self) -> dict[str, Any]
```

Get summary statistics.

**Returns:**

```python
{
    "total_compounds": int,
    "successful": int,
    "failed": int,
    "hazard_class_counts": {"HIGH": n, "MEDIUM": n, "LOW": n},
    "delta_hd_cal_g_min": float | None,
    "delta_hd_cal_g_max": float | None,
}
```

---

## DataFrame Columns

| Column | Type | Description |
|--------|------|-------------|
| `smiles` | str | Input SMILES |
| `canonical_smiles` | str | Canonicalized SMILES |
| `formula` | str | Molecular formula |
| `mw` | float | Molecular weight (g/mol) |
| `delta_hf_kJ_mol` | float | ΔHf° in kJ/mol |
| `delta_hd_kJ_mol` | float | ΔHd in kJ/mol |
| `delta_hd_cal_g` | float | ΔHd in cal/g |
| `ob_percent` | float | Oxygen balance (%) |
| `hazard_class` | str | HIGH, MEDIUM, or LOW |
| `triggered_criteria` | list | CHETAH criteria triggered |
| `alerts` | list | Functional group alerts |
| `gas_volume_L_g` | float | Gas generation (L/g) |
| `error` | str | Error type if failed |
| `error_message` | str | Error details if failed |

---

## Progress Callback

The callback receives `(current_index, total_count)`:

```python
def my_callback(current: int, total: int) -> None:
    print(f"Progress: {current}/{total}")

results = screen(smiles_list, progress_callback=my_callback)
```

### With tqdm

```python
from tqdm import tqdm
from phoenix import screen

pbar = tqdm(total=len(smiles_list))

def update_progress(current, total):
    pbar.update(1)

results = screen(smiles_list, progress_callback=update_progress)
pbar.close()
```

---

## Error Handling

Failed compounds are captured in the DataFrame:

```python
results = screen(smiles_list)
df = results.dataframe

# Find failures
failed = df[df['error'].notna()]
for _, row in failed.iterrows():
    print(f"{row['smiles']}: {row['error']} - {row['error_message']}")
```

### Possible Errors

| Error | Description |
|-------|-------------|
| `InvalidSmilesError` | Invalid SMILES string |
| `UnsupportedElementError` | Unsupported element |
| `UnsupportedStructureError` | Charged/radical species |
| `MissingGroupError` | Missing Benson GA data |
| `UnexpectedError` | Other errors |

---

## Example: Full Workflow

```python
from phoenix import screen

# Read SMILES from file
with open("compounds.smi") as f:
    smiles_list = [line.strip() for line in f if line.strip()]

# Screen with progress
from tqdm import tqdm
pbar = tqdm(total=len(smiles_list), desc="Screening")

results = screen(
    smiles_list,
    progress_callback=lambda c, t: pbar.update(1)
)
pbar.close()

# Get DataFrame
df = results.dataframe

# Summary
summary = results.summary()
print(f"Screened: {summary['total_compounds']}")
print(f"High hazard: {summary['hazard_class_counts'].get('HIGH', 0)}")

# Export
results.to_csv("results.csv")

# Filter
high_risk = df[df['hazard_class'] == 'HIGH']
print(f"\nHigh-risk compounds: {len(high_risk)}")
for _, row in high_risk.iterrows():
    print(f"  {row['formula']}: {row['delta_hd_cal_g']:.0f} cal/g")
```

---

## Large-Scale Processing

### Chunked Processing

```python
import pandas as pd
from phoenix import screen

def screen_in_chunks(smiles_list, chunk_size=1000):
    all_dfs = []
    for i in range(0, len(smiles_list), chunk_size):
        chunk = smiles_list[i:i+chunk_size]
        results = screen(chunk)
        all_dfs.append(results.dataframe)
    return pd.concat(all_dfs, ignore_index=True)
```

### Parallel Processing

```python
from concurrent.futures import ProcessPoolExecutor
from phoenix import screen

def screen_chunk(chunk):
    return screen(chunk).dataframe

def parallel_screen(smiles_list, n_workers=4, chunk_size=100):
    chunks = [smiles_list[i:i+chunk_size]
              for i in range(0, len(smiles_list), chunk_size)]

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        dfs = list(executor.map(screen_chunk, chunks))

    return pd.concat(dfs, ignore_index=True)
```
