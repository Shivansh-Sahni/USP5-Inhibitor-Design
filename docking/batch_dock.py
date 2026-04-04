import subprocess
from pathlib import Path

SMILES_FILE = "SMILES-Top100-04-04-2026.txt"
RECEPTOR = "3IHP.pdbqt"
CONFIG = "config.txt"

OUTDIR = Path("results")
OUTDIR.mkdir(exist_ok=True)

def run(cmd):
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    if result.returncode != 0:
        print("FAILED:", cmd)
        print(result.stderr)
        return False
    return True

def get_score(logfile):
    with open(logfile) as f:
        for line in f:
            parts = line.split()
            if len(parts) > 1 and parts[0].isdigit():
                return parts[1]
    return "NA"

summary = []

with open(SMILES_FILE) as f:
    for i, line in enumerate(f):
        line = line.strip()
        if not line:
            continue

        name = f"cmpd_{i+1}"
        smiles = line

        print("Running", name)

        pdb = OUTDIR / f"{name}.pdb"
        pdb_min = OUTDIR / f"{name}_min.pdb"
        sdf = OUTDIR / f"{name}.sdf"
        pdbqt = OUTDIR / f"{name}.pdbqt"
        out = OUTDIR / f"{name}_out.pdbqt"
        log = OUTDIR / f"{name}.log"

        if not run(f'obabel -:"{smiles}" -O {pdb} --gen3d'):
            continue
        if not run(f'obabel {pdb} -O {pdb_min} --minimize --ff MMFF94'):
            continue
        if not run(f'obabel {pdb_min} -O {sdf}'):
            continue
        if not run(f'mk_prepare_ligand.py -i {sdf} -o {pdbqt}'):
            continue
        if not run(f'vina --receptor {RECEPTOR} --ligand {pdbqt} --config {CONFIG} --out {out} --log {log}'):
            continue

        score = get_score(log)
        summary.append((name, smiles, score))

with open(OUTDIR / "summary.csv", "w") as f:
    f.write("name,smiles,score\n")
    for row in summary:
        f.write(",".join(row) + "\n")

print("DONE")
