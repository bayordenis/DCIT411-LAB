import time
import os
import shutil
import subprocess

# prefer preprocessed file if present, otherwise fall back to original sequences
infile = "data/clean_sequences.fasta"
if not os.path.exists(infile):
    alt = "data/sequences.fasta"
    if os.path.exists(alt):
        print(f"Warning: {infile} not found; falling back to {alt}")
        infile = alt
out_dir = "results/alignments"
os.makedirs(out_dir, exist_ok=True)

results = {}

try:
    from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline, MafftCommandline
    tools = {
        "ClustalW": ClustalwCommandline(input=infile),
        "MUSCLE": MuscleCommandline(input=infile),
        "MAFFT": MafftCommandline(input=infile)
    }

    for name, cline in tools.items():
        start = time.time()
        stdout, stderr = cline()
        end = time.time()
        results[name] = end - start

except Exception:
    # Fallback: try calling common MSA executables directly
    executables = {
        "ClustalW": ["clustalw2", "clustalw"],
        "MUSCLE": ["muscle"],
        "MAFFT": ["mafft"]
    }

    for name, candidates in executables.items():
        exe = None
        for c in candidates:
            if shutil.which(c):
                exe = c
                break
        if not exe:
            print(f"{name} executable not found; skipping.")
            continue

        outfile = os.path.join(out_dir, f"{name.lower()}_aligned.fasta")

        if name == "MUSCLE":
            # modern MUSCLE uses -align / -output flags rather than -in/-out
            cmd = [exe, "-align", infile, "-output", outfile]
            run_stdout = subprocess.DEVNULL
        elif name == "MAFFT":
            # MAFFT writes to stdout by default; redirect into outfile
            cmd = [exe, infile]
            run_stdout = open(outfile, "w")
        else:  # ClustalW
            # clustalw accepts INFILE/OUTFILE as -INFILE=... style on some installs
            cmd = [exe, f"-INFILE={infile}", f"-OUTFILE={outfile}"]
            run_stdout = subprocess.DEVNULL

        start = time.time()
        try:
            if name == "MAFFT":
                subprocess.run(cmd, check=True, stdout=run_stdout, stderr=subprocess.PIPE)
                run_stdout.close()
            else:
                subprocess.run(cmd, check=True, stdout=run_stdout, stderr=subprocess.PIPE)
            end = time.time()
            results[name] = end - start
        except FileNotFoundError:
            print(f"{name} executable not found when running {cmd}; skipping.")
        except subprocess.CalledProcessError as e:
            print(f"{name} failed: returncode={e.returncode}")


print("MSA runtimes (seconds):")
for k, v in results.items():
    print(f"{k}: {v:.3f}")