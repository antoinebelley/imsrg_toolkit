#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import sys
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

# --------------------------
# Input LECs (once)
# --------------------------
LECs = pd.read_csv(
    "/work/submit/abelley/imsrg_toolkit/data/8000Samples.txt",
    usecols=list(range(18))
).set_index("SampleID")

old_LECs_duplicate = [
    40383, 51445, 66697, 97836, 144298, 163571, 172911, 215846, 237612, 252207,
    280489, 28255, 388492, 444839, 537066, 546455, 556707, 557713, 589500, 647468,
    690146, 704538, 708334, 709194, 732716, 777891, 829596
]

# --------------------------
# Column schemas
# --------------------------
LEC_labels = [
    'Ct1S0pp','Ct1S0np','Ct1S0nn','Ct3S1','C1S0','C3P0','C1P1','C3P1','C3S1','CE1','C3P2',
    'c1','c2','c3','c4','cD','cE'
]

# Physics/meta columns kept from the Rp2 row (bra/ket info is the same across rows in examples)
meta_cols = [
    'Nucl bra','J bra','P bra','n bra','Energy bra',
    'Nucl ket','J ket','P ket','n ket','Energy ket'
]

# We keep Zero/One/Two specifically for the Rp2 row (legacy expectations) + new tidy observable columns
obs_scalar_cols = ['Rp2', 'M1', 'E2', 'M1_2BC']  # sums of Zero+One+Two per observable
legacy_cols = ['Zero','One','Two']               # terms for Rp2 specifically
derived_cols = ['Rch']                            # charge radius from Rp2

col_names = meta_cols + legacy_cols + obs_scalar_cols + derived_cols
col_labels = LEC_labels + col_names
ALLOWED_COLS = set(col_labels)  # what we will keep/write, in this exact order

# --------------------------
# Periodic table (for Z lookup)
# --------------------------
PERIODIC = {
    1:"H",2:"He",3:"Li",4:"Be",5:"B",6:"C",7:"N",8:"O",9:"F",10:"Ne",
    11:"Na",12:"Mg",13:"Al",14:"Si",15:"P",16:"S",17:"Cl",18:"Ar",19:"K",20:"Ca",
    21:"Sc",22:"Ti",23:"V",24:"Cr",25:"Mn",26:"Fe",27:"Co",28:"Ni",29:"Cu",30:"Zn",
    31:"Ga",32:"Ge",33:"As",34:"Se",35:"Br",36:"Kr",37:"Rb",38:"Sr",39:"Y",40:"Zr",
    41:"Nb",42:"Mo",43:"Tc",44:"Ru",45:"Rh",46:"Pd",47:"Ag",48:"Cd",49:"In",50:"Sn",
    51:"Sb",52:"Te",53:"I",54:"Xe",55:"Cs",56:"Ba",57:"La",58:"Ce",59:"Pr",60:"Nd",
    61:"Pm",62:"Sm",63:"Eu",64:"Gd",65:"Tb",66:"Dy",67:"Ho",68:"Er",69:"Tm",70:"Yb",
    71:"Lu",72:"Hf",73:"Ta",74:"W",75:"Re",76:"Os",77:"Ir",78:"Pt",79:"Au",80:"Hg",
    81:"Tl",82:"Pb",83:"Bi",84:"Po",85:"At",86:"Rn",87:"Fr",88:"Ra",89:"Ac",90:"Th",
    91:"Pa",92:"U",93:"Np",94:"Pu",95:"Am",96:"Cm",97:"Bk",98:"Cf",99:"Es",100:"Fm",
    101:"Md",102:"No",103:"Lr",104:"Rf",105:"Db",106:"Sg",107:"Bh",108:"Hs",109:"Mt",
    110:"Ds",111:"Rg",112:"Cn",113:"Nh",114:"Fl",115:"Mc",116:"Lv",117:"Ts",118:"Og"
}
SYM2Z = {v: k for k, v in PERIODIC.items()}

# --------------------------
# Helpers
# --------------------------
def get_file_folder_list(base_dir: str, pattern: str):
    """Return subfolders under base_dir matching pattern like 'Ca*'."""
    folder_list = glob.glob(os.path.join(base_dir, pattern))
    folder_list = [f for f in folder_list if os.path.isdir(f)]
    return folder_list

def Rp2_to_Rch2(Rp2, Z, N, CODATA=True):
    """
    Convert mean squared point proton radius to mean squared charge radius.
    """
    if CODATA:
        rcp2 = 0.8783**2  # CODATA
        rcn2 = -0.115     # CODATA
    else:
        rcp2 = 0.709      # Nature 466, 213 (2010)
        rcn2 = -0.106     # Phys. Rev. Lett. 124, 082501

    DF = 0.033
    return Rp2 + rcp2 + N/Z * rcn2 + DF

def read_0b_energy(fn_csv):
    """Read 0b term from paired .snt file (unused here but kept)."""
    fn = fn_csv.replace(".csv", ".snt")
    with open(fn, 'r') as f:
        lines = f.readlines()
    import re as _re
    zerob = float(_re.findall(r'[+-]?[0-9]+.+[0-9]', lines[4])[0])
    return zerob

# --------------------------
# Results dataframe wrapper
# --------------------------
class IMSRGResultsDF:
    """
    Data frame object containing all the results up to date for the decay.
    Index: (Sample, emax)
    Columns: LEC_labels + physics result columns
    """
    def __init__(self, file=None):
        self.intialize_df(file)

    def create_df(self):
        names = ['Sample', 'emax']
        empty_index = pd.MultiIndex.from_tuples([], names=names)
        self.df = pd.DataFrame(index=empty_index, columns=col_labels)

    def intialize_df(self, file=None):
        if file is None:
            self.create_df()
        else:
            try:
                df = pd.read_csv(file, dtype={'Sample': int, 'emax': int})
                df.set_index(['Sample', 'emax'], inplace=True)
                df.sort_index(inplace=True)
                # ensure correct column order (and drop any unexpected columns)
                df = df.reindex(columns=col_labels)
                self.df = df
                print("Read df from file.")
            except Exception:
                print("Output csv file not found or unreadable. Creating new dataframe.")
                self.create_df()

    def add_to_dataframe(self, Nucl, Sample, emax, data: dict):
        index = (Sample, emax)

        # filter to allowed columns only (drop 'fn_op' or any stray fields)
        filtered = {k: v for k, v in data.items() if k in ALLOWED_COLS}

        # normalize strings
        for k, v in list(filtered.items()):
            if isinstance(v, str):
                filtered[k] = v.strip()

        # write values
        for key, val in filtered.items():
            try:
                if self.df.loc[index, key] == val:
                    continue
            except KeyError:
                pass
            self.df.loc[index, key] = val

        # housekeeping + enforce canonical col order
        self.df = self.df.fillna('')
        self.df.drop_duplicates(inplace=True)
        self.df.sort_index(inplace=True)
        self.df = self.df.reindex(columns=col_labels)

    def to_csv(self, file):
        self.df.sort_index(inplace=True)
        self.df.drop_duplicates(inplace=True)
        # enforce final column order just before writing
        self.df = self.df.reindex(columns=col_labels)
        self.df.to_csv(file)

    def __str__(self):
        return self.df.__str__()

# --------------------------
# CSV parsing helpers
# --------------------------
_OBS_PATTERN = re.compile(r'_(Rp2|M1|E2|M1_2BC)\.snt$')

def _read_result_rows(csv_path: str) -> pd.DataFrame:
    """
    Read the CSV robustly, keeping 'fn_op' so we can identify the observable.
    Handles the legacy single-row case and multi-row case.
    """
    # Load first, then normalize columns and drop an unnamed index column if present
    df = pd.read_csv(csv_path)
    # normalize headers
    df.columns = [str(c).strip() for c in df.columns]
    # drop pandas' default index column if present
    if df.columns[0].lower().startswith('unnamed'):
        df = df.drop(columns=[df.columns[0]])

    # in some legacy outputs, 'fn_op' may be missing; we still can handle one-row Rp2
    return df

def _extract_observable_from_fn(fn_op: str) -> str:
    """
    From the 'fn_op' path, extract the observable label among ['Rp2','M1','E2','M1_2BC'].
    Returns '' if not found.
    """
    if not isinstance(fn_op, str):
        return ''
    m = _OBS_PATTERN.search(fn_op)
    return m.group(1) if m else ''

# --------------------------
# Core processing
# --------------------------
def process_nuclei(Nucl, Z, N, directory_glob, output_file):
    """
    Process all CSVs under directory_glob (e.g., '/path/Ca49/*/*.csv') for a given isotope.
    """
    results = IMSRGResultsDF(file=output_file)

    for name in tqdm(glob.glob(os.path.join(directory_glob, "*.csv"))):
        fn = os.path.basename(name)
        params = fn.split("_")

        SampleId, emax = None, None
        for param in params:
            if param.isnumeric():
                SampleId = int(param)
            elif param.startswith('e'):
                try:
                    emax = int(param.replace('e', ''))
                except ValueError:
                    pass

        if SampleId is None or emax is None:
            # skip files that don't match the expected naming convention
            continue

        if SampleId in old_LECs_duplicate:
            continue

        # Pull LECs for this sample
        try:
            LEC_cols = LECs.loc[SampleId].to_dict()
        except KeyError:
            # If SampleId not in LECs table, skip gracefully
            continue

        # --------------------------
        # Read physics rows
        # --------------------------
        try:
            df_rows = _read_result_rows(name)
        except Exception:
            # Fallback to the older, stricter reader (rarely needed now)
            try:
                dict_cols = pd.read_csv(
                    name,
                    dtype={
                        'Nucl bra': str, 'J bra': str, 'P bra': str, 'n bra': str, 'Energy bra': float,
                        'Nucl ket': str, 'J ket': str, 'P ket': str, 'n ket': str, 'Energy ket': float,
                        'Zero': float, 'One': float, 'Two': float
                    },
                    usecols=[i for i in range(16)][2:]  # skip col 0 (df index) & col 1 (fn_op)
                )
                df_rows = dict_cols
            except Exception:
                # Give up on this file
                continue

        # Normalize column names and types
        df_rows.columns = [c.strip() for c in df_rows.columns]
        # Ensure required numeric fields are present when they appear
        for num_col in ['Energy bra','Energy ket','Zero','One','Two']:
            if num_col in df_rows.columns:
                df_rows[num_col] = pd.to_numeric(df_rows[num_col], errors='coerce')

        # --------------------------
        # Aggregate per observable
        # --------------------------
        # Initialize output payload with zeros for the scalar obs columns
        out = {
            'Rp2': 0.0, 'M1': 0.0, 'E2': 0.0, 'M1_2BC': 0.0,
            'Zero': '', 'One': '', 'Two': '',  # legacy (we keep these for Rp2 only)
            'Rch': ''
        }

        # We also want to keep bra/ket/meta columns (take from the first row)
        for k in meta_cols:
            out[k] = ''

        # Determine if we are in legacy single-row mode or multi-row mode
        rows_count = len(df_rows)

        if rows_count == 0:
            # nothing to do
            continue

        if rows_count == 1:
            # Legacy: single Rp2 entry, possibly without 'fn_op'
            row = df_rows.iloc[0].to_dict()
            # Fill meta (if present)
            for k in meta_cols:
                if k in row:
                    out[k] = row[k]

            # Compute Rp2 from terms
            try:
                zero = float(row.get('Zero', 0.0))
                one  = float(row.get('One',  0.0))
                two  = float(row.get('Two',  0.0))
            except Exception:
                zero, one, two = 0.0, 0.0, 0.0

            Rp2_val = zero + one + two
            out['Rp2'] = Rp2_val
            out['Zero'] = zero
            out['One']  = one
            out['Two']  = two
            # Others default 0.0 (already set)
            # Compute Rch
            out['Rch'] = np.sqrt(Rp2_to_Rch2(Rp2_val, Z, N))

        else:
            # Multi-row: parse each row's observable from 'fn_op'
            # Take meta from the first row (they should be identical across rows for the same file)
            first = df_rows.iloc[0].to_dict()
            for k in meta_cols:
                if k in first:
                    out[k] = first[k]

            # Iterate rows and fill observable sums
            seen_rp2_terms = False
            for _, r in df_rows.iterrows():
                rdict = r.to_dict()
                obs = _extract_observable_from_fn(str(rdict.get('fn_op', '')))
                # If fn_op is missing or unparsable, try best-effort by row order
                if obs == '':
                    # Assume first row is Rp2, then M1, E2, M1_2BC by order if available
                    # Map by index position
                    idx = _
                    order = ['Rp2', 'M1', 'E2', 'M1_2BC']
                    if isinstance(idx, (int, np.integer)) and 0 <= idx < len(order):
                        obs = order[int(idx)]
                    else:
                        continue  # skip if we truly can't tell

                try:
                    zero = float(rdict.get('Zero', 0.0))
                    one  = float(rdict.get('One',  0.0))
                    two  = float(rdict.get('Two',  0.0))
                except Exception:
                    zero, one, two = 0.0, 0.0, 0.0

                val = zero + one + two

                if obs == 'Rp2':
                    out['Rp2'] = val
                    # preserve term breakdown specifically for Rp2 (legacy columns)
                    out['Zero'] = zero
                    out['One']  = one
                    out['Two']  = two
                    seen_rp2_terms = True
                elif obs in out:
                    out[obs] = val

            # If we never saw an Rp2 row, leave Rch blank; otherwise compute Rch
            if seen_rp2_terms:
                out['Rch'] = np.sqrt(Rp2_to_Rch2(out['Rp2'], Z, N))
            else:
                # No Rp2 present; keep Rch empty string (or set to NaN)
                out['Rch'] = ''

        # Merge LECs and filtered physics/meta into one record
        record = {}
        record.update(LEC_cols)
        record.update(out)

        # Final filter to allowed columns only
        record = {k: v for k, v in record.items() if k in ALLOWED_COLS}

        # Store
        results.add_to_dataframe(Nucl, SampleId, emax, record)

    return results

def folder_name_to_info(folder_path):
    """
    Input:  /home/submit/josemm/results/Ca49
    Output: ('Ca', Z=20, N=A-Z) with A=49
    """
    # Remove 'partial_' prefix if present
    folder_path = re.sub(r'partial_', '', folder_path)
    tail = os.path.basename(folder_path.rstrip("/"))
    sym = re.sub(r'[^A-Za-z]', '', tail)  # 'Ca'
    A_str = re.sub(r'[^0-9]', '', tail)   # '49'
    if not A_str:
        raise ValueError(f"Cannot parse mass number from: {folder_path}")
    A = int(A_str)
    if sym not in SYM2Z:
        raise ValueError(f"Unknown element symbol '{sym}' in: {folder_path}")
    z = SYM2Z[sym]
    n = A - z
    return sym, z, n, A

def get_output_file(out_dir, sym, A):
    nuclide = f"{sym}{A}"
    return os.path.join(out_dir, f"{nuclide}_results.csv")

# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    PATH_2_RES = "/home/submit/josemm/results"
    # PATH_2_RES = "/ceph/submit/data/user/j/josemm/IMSRG_RES/results/CALCIUM"
    PATH_4_OUTPUT = "/work/submit/josemm/imsrg_results"
    PATTERN = "Ca*"   # e.g., 'Ca*' to process all Calcium isotopes

    # find top-level isotope folders like /results/Ca49, /results/Ca48, ...
    folders = get_file_folder_list(PATH_2_RES, PATTERN)
    print(folders)

    for folder in folders:
        if "Ca47" not in folder:
            continue  # TEMP: process only partials for testing
        print(f"Found folder: {folder}")
        try:
            sym, Z, N, A = folder_name_to_info(folder)
            print(f"Processing isotope: {sym}{A} (Z={Z}, N={N})")
        except ValueError as e:
            print(f"Skipping '{folder}': {e}")
            continue

        directory_glob = os.path.join(folder, "**")
        output_file = get_output_file(PATH_4_OUTPUT, sym, A)

        results = process_nuclei(sym, Z, N, directory_glob, output_file)
        results.to_csv(output_file)
        print(f"Processed {sym}{A} (Z={Z}, N={N}) and saved to {output_file}")
    print("All done.")
