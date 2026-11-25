
import pandas as pd
import numpy as np
import os
import re

# =====================================================
# 1. File path
# =====================================================

from pathlib import Path
import argparse
_root = Path(__file__).resolve().parent
_parser = argparse.ArgumentParser(description='Convert the supplied (y,x,t) Excel input to NPZ.')
_parser.add_argument('--input', type=Path, default=_root/'data/2D_RS_Real_1_t.xlsx')
_parser.add_argument('--output', type=Path, default=_root/'outputs/2D_RS_Real_1_t.npz')
_args = _parser.parse_args()
file_path = str(_args.input.resolve())

if not os.path.exists(file_path):
    raise FileNotFoundError(
        f"Cannot find file:\n{file_path}"
    )

# =====================================================
# 2. Read Excel workbook
# =====================================================

xl = pd.ExcelFile(file_path)
sheet_names = xl.sheet_names

print(f"Number of sheets: {len(sheet_names)}")
print(f"Sheet names: {sheet_names}")


# =====================================================
# 3. Extract t value from sheet name
# =====================================================

def parse_t(sheet_name):

    m = re.search(r'[-+]?\d*\.?\d+', sheet_name)

    if m is None:
        raise ValueError(
            f"Cannot find t value in sheet name: {sheet_name}"
        )

    return float(m.group())


# =====================================================
# 4. Read data from each sheet
# =====================================================

records = []

for name in sheet_names:

    print(f"\nProcessing sheet: {name}")

    t_val = parse_t(name)

    # IMPORTANT:
    # First row = x coordinates
    # First column = y coordinates
    raw_df = pd.read_excel(
        xl,
        sheet_name=name,
        header=0,
        index_col=0
    )

    # Convert column names and row indices to floats
    x_values = np.round(
        pd.to_numeric(raw_df.columns).to_numpy(dtype=float),
        6
    )

    y_values = np.round(
        pd.to_numeric(raw_df.index).to_numpy(dtype=float),
        6
    )

    # Check data
    values = raw_df.to_numpy(dtype=float)

    if not np.all(np.isfinite(values)):
        raise ValueError(
            f"Missing or non-finite values in sheet: {name}"
        )

    # Construct coordinate grids
    x_grid, y_grid = np.meshgrid(
        x_values,
        y_values
    )

    # Convert to long format
    long_df = pd.DataFrame({
        'y': y_grid.flatten(),
        'x': x_grid.flatten(),
        'value': values.flatten()
    })

    long_df['t'] = t_val

    records.append(long_df)


# =====================================================
# 5. Combine all time slices
# =====================================================

long_df = pd.concat(
    records,
    ignore_index=True
)


# =====================================================
# 6. Extract unique coordinates and time
# =====================================================

t_arr = np.sort(
    long_df['t'].unique()
)

x_arr = np.sort(
    long_df['x'].unique()
)

y_arr = np.sort(
    long_df['y'].unique()
)

print("\n========== Coordinate Information ==========")

print(f"t coordinates: {t_arr}")
print(f"x coordinates: {x_arr}")
print(f"y coordinates: {y_arr}")

print(f"Number of t: {len(t_arr)}")
print(f"Number of x: {len(x_arr)}")
print(f"Number of y: {len(y_arr)}")


# =====================================================
# 7. Construct solution array
# =====================================================

# Shape: (number of y, number of x, number of t)

u_sol = np.zeros(
    (len(y_arr), len(x_arr), len(t_arr))
)

for i, t in enumerate(t_arr):

    slice_df = long_df[
        long_df['t'] == t
    ]

    pivot_df = slice_df.pivot(
        index='y',
        columns='x',
        values='value'
    )

    # Ensure consistent coordinate ordering
    pivot_df = pivot_df.reindex(
        index=y_arr,
        columns=x_arr
    )

    if pivot_df.isna().any().any():
        raise ValueError(
            f"Missing grid points at t = {t}"
        )

    u_sol[:, :, i] = pivot_df.to_numpy()


# =====================================================
# 8. Save as NPZ file
# =====================================================

save_path = str(_args.output.resolve())
Path(save_path).parent.mkdir(parents=True, exist_ok=True)

np.savez(
    save_path,
    t=t_arr.reshape(-1, 1),
    x=x_arr.reshape(-1, 1),
    y=y_arr.reshape(-1, 1),
    u_sol=u_sol
)


# =====================================================
# 9. Print final information
# =====================================================

print("\n========== Final Results ==========")

print(f"t array:\n{t_arr}")

print(f"\nx array:\n{x_arr}")

print(f"\ny array:\n{y_arr}")

print(f"\nu_sol shape: {u_sol.shape} (y, x, t)")

print(f"\nFile saved successfully:\n{save_path}")

print("\nAll done!")
