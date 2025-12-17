#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
# AlgPred 3.0 is developed for predicting allergenic and non-allergenic      #
# proteins/peptides from their primary sequence. It is developed by Prof     #
# G. P. S. Raghava's group. Please cite: AlgPred 3.0                          #
##############################################################################

import argparse
import os
import sys
import warnings
from pathlib import Path

import pandas as pd
import joblib

warnings.filterwarnings("ignore")

# ===================== PATHS & CONSTANTS =====================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MODEL_PATH = "algpred3_model.sav"

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
VALID_AAS = set(AA_LIST)
LOG_FILE = "stand_error.log"

# ===================== ARGUMENT PARSER =====================
def parse_args():
    parser = argparse.ArgumentParser(
        description="AlgPred 3.0 unified script: prediction (pred), protein scan (scan), design (des)"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA or plain text file")
    parser.add_argument("-o", "--output", default="final_predictions.csv",
                        help="Output CSV file")
    parser.add_argument("-j", "--job", required=True,
                        choices=["pred", "scan", "des"],
                        help="Job type")
    parser.add_argument("-t", "--threshold", type=float, default=0.5,
                        help="Prediction threshold")

    # Protein Scan specific
    parser.add_argument("-l", "--length", type=int, default=None,
                        help="Sliding window length (Protein Scan mode)")
    parser.add_argument("-s", "--step", type=int, default=1,
                        help="Sliding window step size (Protein Scan mode)")

    parser.add_argument("-wd", "--working_directory", default=None,
                        help="Working directory")

    return parser.parse_args()

# ===================== FASTA UTILITIES =====================
def is_fasta_file(path):
    with open(path) as f:
        for line in f:
            if line.strip():
                return line.startswith(">")
    return False


def read_fasta(path):
    records = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if header:
            records.append((header, "".join(seq)))
    return records


def write_fasta(records, path):
    with open(path, "w") as f:
        for h, s in records:
            f.write(f">{h}\n{s}\n")

# ===================== CLEAN INPUT =====================
def clean_fasta(input_path):
    cleaned = f"{Path(input_path).stem}_clean.fasta"
    removed = []
    kept = []

    fasta_mode = is_fasta_file(input_path)
    print(f"Input detected as {'FASTA' if fasta_mode else 'plain text'}")

    if fasta_mode:
        records = read_fasta(input_path)
    else:
        records = [
            (f"seq{i+1}", line.strip())
            for i, line in enumerate(open(input_path))
            if line.strip()
        ]

    for name, seq in records:
        seq = seq.upper().replace(" ", "")
        bad = sorted(set(seq) - VALID_AAS)
        if bad or not seq:
            removed.append((name, bad if bad else ["<empty>"]))
        else:
            kept.append((name, seq))

    if removed:
        print("⚠️ Removed sequences:")
        with open(LOG_FILE, "a") as log:
            log.write("\n===== REMOVED SEQUENCES =====\n")
            for name, bad in removed:
                msg = f"{name}: {', '.join(bad)}"
                print("  " + msg)
                log.write(msg + "\n")

    if not kept:
        print("❌ No valid sequences remain.")
        return None

    write_fasta(kept, cleaned)
    print(f"✅ Cleaned FASTA saved → {cleaned} ({len(kept)} valid sequences)")
    return cleaned

# ===================== FEATURE EXTRACTION =====================
def dpc_for_sequence(seq):
    dpc = {f"DPC1_{a}{b}": 0 for a in AA_LIST for b in AA_LIST}
    total = len(seq) - 1
    if total <= 0:
        return dpc
    for i in range(total):
        dpc[f"DPC1_{seq[i]}{seq[i+1]}"] += 1
    for k in dpc:
        dpc[k] /= total
    return dpc


def extract_dpc(fasta, output_csv):
    records = read_fasta(fasta)
    df = pd.DataFrame([dpc_for_sequence(seq) for _, seq in records])
    df.to_csv(output_csv, index=False)
    print(f"✅ DPC features saved → {output_csv} ({len(df)} sequences)")
    return df

# ===================== MODEL =====================
def load_model():
    if not os.path.exists(DEFAULT_MODEL_PATH):
        sys.exit(f"❌ Model file missing: {DEFAULT_MODEL_PATH}")
    print(f"Loading model: {DEFAULT_MODEL_PATH}")
    return joblib.load(DEFAULT_MODEL_PATH)


def predict(model, feature_csv, threshold):
    X = pd.read_csv(feature_csv)
    probs = model.predict_proba(X)[:, 1]
    preds = (probs >= threshold).astype(int)
    return probs, preds

# ===================== PROTEIN SCAN UTILITIES =====================
def sliding_windows(name, seq, win_len, step):
    windows = []
    if len(seq) < win_len:
        return windows

    for i in range(0, len(seq) - win_len + 1, step):
        windows.append({
            "ParentSeq": name,
            "Start": i + 1,
            "End": i + win_len,
            "Peptide": seq[i:i + win_len]
        })
    return windows

# ===================== JOB MODES =====================
def run_scan(cleaned_fasta, model, threshold, length, step, output):
    print("=== PROTEIN SCAN MODE (Sliding Window Prediction) ===")

    if length is None:
        sys.exit("❌ Protein Scan mode requires --length")

    records = read_fasta(cleaned_fasta)
    windows = []

    for name, seq in records:
        windows.extend(sliding_windows(name, seq, length, step))

    if not windows:
        sys.exit("❌ No windows generated")

    win_fasta = "scan_windows.fasta"
    write_fasta(
        [(f"{w['ParentSeq']}_{w['Start']}_{w['End']}", w["Peptide"]) for w in windows],
        win_fasta
    )

    feature_csv = "scan_windows_DPC.csv"
    extract_dpc(win_fasta, feature_csv)

    probs, preds = predict(model, feature_csv, threshold)

    df = pd.DataFrame({
        "ParentSeq": [w["ParentSeq"] for w in windows],
        "Start":     [w["Start"] for w in windows],
        "End":       [w["End"] for w in windows],
        "Peptide":   [w["Peptide"] for w in windows],
        "Score":     probs,
        "Prediction": ["Allergen" if p else "Non-Allergen" for p in preds]
    })

    df.to_csv(output, index=False)
    print(f"✅ Protein Scan results saved → {output}")


def run_pred(cleaned_fasta, model, threshold, output):
    print("=== PRED MODE ===")
    feature_csv = f"{Path(cleaned_fasta).stem}_DPC.csv"
    extract_dpc(cleaned_fasta, feature_csv)

    probs, preds = predict(model, feature_csv, threshold)
    names = [h for h, _ in read_fasta(cleaned_fasta)]

    df = pd.DataFrame({
        "Sequence_ID": names,
        "Probability": probs,
        "Status": ["Allergen" if p else "Non-Allergen" for p in preds]
    })
    df.to_csv(output, index=False)
    print(f"✅ Predictions saved → {output}")


def run_des(cleaned_fasta, model, threshold, output):
    print("=== DES MODE ===")

    mutants = []
    originals = read_fasta(cleaned_fasta)

    for name, seq in originals:
        mutants.append((name, "original", seq))
        for i, aa in enumerate(seq):
            for m in AA_LIST:
                if m != aa:
                    mutants.append((name, f"mut{i+1}{aa}>{m}",
                                    seq[:i] + m + seq[i+1:]))

    mut_fasta = f"{Path(cleaned_fasta).stem}_mutants.fasta"
    write_fasta([(f"{sid}_{mid}", s) for sid, mid, s in mutants], mut_fasta)

    feature_csv = f"{Path(mut_fasta).stem}_DPC.csv"
    extract_dpc(mut_fasta, feature_csv)

    probs, preds = predict(model, feature_csv, threshold)

    df = pd.DataFrame({
        "SeqID": [sid for sid, _, _ in mutants],
        "MutantID": [mid for _, mid, _ in mutants],
        "Sequence": [s for _, _, s in mutants],
        "Score": probs,
        "Prediction": ["Allergen" if p else "Non-Allergen" for p in preds]
    })

    df.to_csv(output, index=False)
    print(f"✅ Design results saved → {output}")

# ===================== MAIN =====================
def main():
    print("##############################################################################")
    print("# AlgPred 3.0                                                                #")
    print("#                                                                            #")
    print("# Prediction of Allergenic and Non-Allergenic Proteins/Peptides             #")
    print("#                                                                            #")
    print("# Developed by:                                                              #")
    print("# Raghava Group                                                              #")
    print("# Department of Computational Biology                                        #")
    print("# Indraprastha Institute of Information Technology, Delhi (IIIT-Delhi)       #")
    print("#                                                                            #")
    print("# Please cite:                                                               #")
    print("# AlgPred 3.0                                                                #")
    print("#                                                                            #")
    print("##############################################################################\n")


    args = parse_args()
    os.chdir(args.working_directory or SCRIPT_DIR)

    open(LOG_FILE, "w").close()

    cleaned = clean_fasta(args.input)
    if not cleaned:
        sys.exit(1)

    model = load_model()

    if args.job == "scan":
        run_scan(cleaned, model, args.threshold,
                 args.length, args.step, args.output)
    elif args.job == "pred":
        run_pred(cleaned, model, args.threshold, args.output)
    elif args.job == "des":
        run_des(cleaned, model, args.threshold, args.output)

    print("\n======= Thanks for using AlgPred 3.0 =======\n")

# ===================== ENTRY POINT =====================
if __name__ == "__main__":
    main()
