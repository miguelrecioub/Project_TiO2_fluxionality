import sys
import os
import numpy as np
import pandas as pd

def read_xyz_trajectory(filename):
    """Read multi-frame XYZ; return (elements, coords[T,F,3])."""
    elements = None
    frames = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    i = 0
    n = len(lines)
    while i < n:
        try:
            n_atoms = int(lines[i].strip())
        except Exception:
            break
        block = lines[i + 2 : i + 2 + n_atoms]
        if len(block) < n_atoms:
            break
        elems, xyz = [], []
        for ln in block:
            parts = ln.split()
            if len(parts) < 4:
                raise ValueError(f"Malformed XYZ line: {ln}")
            elems.append(parts[0])
            xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])
        if elements is None:
            elements = elems
        elif elems != elements:
            raise ValueError("Atom ordering or elements change between frames.")
        frames.append(np.array(xyz, dtype=float))
        i += n_atoms + 2
    coords = np.stack(frames, axis=0)
    return elements, coords

def compute_rmsf(coords):
    """Compute per-atom RMSF from pre-aligned coordinates."""
    mean_pos = coords.mean(axis=0)
    disp = coords - mean_pos
    sq = (disp ** 2).sum(axis=2)
    return np.sqrt(sq.mean(axis=0))

def compute_sigma_per_element(elements, rmsf):
    """Compute σ_X and RMSF spread for each element."""
    df = pd.DataFrame({"Element": elements, "RMSF(Å)": rmsf})
    results = []
    for elem, sub in df.groupby("Element"):
        sigma = np.sqrt(np.mean(sub["RMSF(Å)"].values ** 2))
        spread = np.std(sub["RMSF(Å)"].values, ddof=1) if len(sub) > 1 else 0.0
        results.append({
            "Element": elem,
            "σ_X (Å)": sigma,
            "Std_RMSF_across_atoms (Å)": spread,
            "N_atoms": len(sub)
        })
    return pd.DataFrame(results)

def main():
    if len(sys.argv) != 2:
        print("Usage: python compute_sigma_per_element.py trajectory.xyz")
        sys.exit(1)
    xyz_file = sys.argv[1]
    if not os.path.exists(xyz_file):
        print(f"Error: file '{xyz_file}' not found.")
        sys.exit(1)

    print(f"\n Reading trajectory: {xyz_file}")
    elements, coords = read_xyz_trajectory(xyz_file)
    rmsf = compute_rmsf(coords)
    df = compute_sigma_per_element(elements, rmsf)

    out_csv = os.path.splitext(xyz_file)[0] + "_sigma_per_element.csv"
    df.to_csv(out_csv, index=False)
    print("\n RMSF per element (σ_X) results:")
    print(df.to_string(index=False))
    print(f"\n Results saved to: {out_csv}\n")

if __name__ == "__main__":
    main()
