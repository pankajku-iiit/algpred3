#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##############################################################################
# AlgPred 3.0
# Prediction of Allergenic and Non-Allergenic Peptides
# Developed by Raghava Group, IIIT-Delhi
##############################################################################

import argparse
import os
import sys
import warnings
import subprocess
from pathlib import Path
import urllib.request
import pandas as pd
import joblib

warnings.filterwarnings("ignore")

# ===================== PATHS & CONSTANTS =====================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

DEFAULT_MODEL_PATH = os.path.join(SCRIPT_DIR, "algpred3_model.sav")
PFEATURE_BIN = os.path.join(SCRIPT_DIR, "pfeature_comp")
COLUMNS_FILE = os.path.join(SCRIPT_DIR, "columns.csv")

MODEL_URL = "https://github.com/pankajku-iiit/algpred3/releases/download/v3.0.0/algpred3_model.sav"
PFEATURE_URL = "https://github.com/pankajku-iiit/algpred3/releases/download/v3.0.0/pfeature_comp"
COLUMNS_URL = "https://github.com/pankajku-iiit/algpred3/releases/download/v3.0.0/columns.csv"

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
VALID_AAS = set(AA_LIST)
LOG_FILE = "stand_error.log"
MODEL_URL

# ===================== ARGUMENT PARSER =====================
def parse_args():
    parser = argparse.ArgumentParser(
        description="AlgPred 3.0 unified script: pred / scan / des"
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA or plain text file")
    parser.add_argument("-o", "--output", default="final_predictions.csv",
                        help="Output CSV")
    parser.add_argument("-j", "--job", required=True,
                        choices=["pred", "scan", "des"])
    parser.add_argument("-t", "--threshold", type=float, default=0.5)

    parser.add_argument("-l", "--length", type=int,
                        help="Window length (scan mode)")
    parser.add_argument("-s", "--step", type=int, default=1)

    return parser.parse_args()


# ===================== AUTO DOWNLOAD =====================
def download_file(url, path, make_executable=False):
    print(f"⬇ Downloading {os.path.basename(path)}...")
    try:
        urllib.request.urlretrieve(url, path)
        if make_executable:
            os.chmod(path, 0o755)
        print(f"✅ {os.path.basename(path)} ready")
    except Exception as e:
        sys.exit(f"❌ Failed to download {path}: {e}")


def setup_dependencies():
    if not os.path.exists(PFEATURE_BIN):
        download_file(PFEATURE_URL, PFEATURE_BIN, make_executable=True)

    if not os.path.exists(COLUMNS_FILE):
        download_file(COLUMNS_URL, COLUMNS_FILE)

    if not os.path.exists(DEFAULT_MODEL_PATH):
        download_file(MODEL_URL, DEFAULT_MODEL_PATH)


# ===================== FASTA UTILITIES =====================
def is_fasta_file(path):
    with open(path) as f:
        for line in f:
            if line.strip():
                return line.startswith(">")
    return False


def read_fasta(path):
    records, header, seq = [], None, []

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
    print("\n🔍 Checking input sequences...")

    cleaned = f"{Path(input_path).stem}_clean.fasta"
    removed, kept = [], []

    fasta_mode = is_fasta_file(input_path)
    print(f"Input detected as {'FASTA' if fasta_mode else 'Plain Text'}")

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
        print("⚠ Removed sequences:")
        with open(LOG_FILE, "a") as log:
            log.write("\n===== REMOVED SEQUENCES =====\n")
            for name, bad in removed:
                msg = f"{name}: {','.join(bad)}"
                print("  " + msg)
                log.write(msg + "\n")

    if not kept:
        print("❌ No valid sequences remain.")
        return None

    write_fasta(kept, cleaned)
    print(f"✅ {len(kept)} valid sequences saved → {cleaned}")

    return cleaned


# ===================== FEATURE EXTRACTION =====================
def extract_features(input_fasta):
    print("⚙ Running pfeature...")

    temp_csv = "temp_features.csv"

    cmd = [
        PFEATURE_BIN,
        "-i", input_fasta,
        "-j", "AAC,DPC,BTC,DDR,PAAC,APAAC,CeTD",
        "-o", temp_csv
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"❌ pfeature failed: {e}")

    print("🔬 Selecting trained feature columns...")

    try:
        df = pd.read_csv(temp_csv)
        required_cols = pd.read_csv(COLUMNS_FILE)["Columns"].tolist()

        missing = set(required_cols) - set(df.columns)
        if missing:
            sys.exit(f"❌ Missing columns: {missing}")

        df_selected = df[required_cols]

    except Exception as e:
        sys.exit(f"❌ Feature processing failed: {e}")

    os.remove(temp_csv)

    print(f"✅ Feature matrix ready ({df_selected.shape[0]} × {df_selected.shape[1]})")

    return df_selected


# ===================== MODEL =====================
def load_model():
    print("📦 Loading trained model...")
    return joblib.load(DEFAULT_MODEL_PATH)


def predict(model, X_df, threshold):
    print("🤖 Running prediction...")
    probs = model.predict_proba(X_df)[:, 1]
    preds = (probs >= threshold).astype(int)
    return probs, preds


# ===================== PROTEIN SCAN =====================
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


def run_scan(cleaned_fasta, model, threshold, length, step, output):
    print("\n=== PROTEIN SCAN MODE ===")

    if length is None:
        sys.exit("❌ Provide window length for scan mode")

    records = read_fasta(cleaned_fasta)
    windows = []

    for name, seq in records:
        windows.extend(sliding_windows(name, seq, length, step))

    print(f"🧬 Generated {len(windows)} sliding windows")

    win_fasta = "scan_windows.fasta"
    write_fasta(
        [(f"{w['ParentSeq']}_{w['Start']}_{w['End']}", w["Peptide"])
         for w in windows],
        win_fasta
    )

    X = extract_features(win_fasta)
    probs, preds = predict(model, X, threshold)

    df = pd.DataFrame({
        "ParentSeq": [w["ParentSeq"] for w in windows],
        "Start": [w["Start"] for w in windows],
        "End": [w["End"] for w in windows],
        "Peptide": [w["Peptide"] for w in windows],
        "Score": probs,
        "Prediction": ["Allergen" if p else "Non-Allergen" for p in preds]
    })

    df.to_csv(output, index=False)
    print(f"✅ Scan results saved → {output}")


# ===================== PRED MODE =====================
def run_pred(cleaned_fasta, model, threshold, output):
    print("\n=== PRED MODE ===")

    X = extract_features(cleaned_fasta)
    probs, preds = predict(model, X, threshold)

    names = [h for h, _ in read_fasta(cleaned_fasta)]

    df = pd.DataFrame({
        "Sequence_ID": names,
        "Probability": probs,
        "Status": ["Allergen" if p else "Non-Allergen" for p in preds]
    })

    df.to_csv(output, index=False)
    print(f"✅ Predictions saved → {output}")


# ===================== DESIGN MODE =====================
def run_des(cleaned_fasta, model, threshold, output):
    print("\n=== DESIGN MODE ===")

    mutants = []
    originals = read_fasta(cleaned_fasta)

    for name, seq in originals:
        mutants.append((name, "original", seq))
        for i, aa in enumerate(seq):
            for m in AA_LIST:
                if m != aa:
                    mutants.append((name, f"mut{i+1}{aa}-{m}",
                                    seq[:i] + m + seq[i+1:]))

    print(f"🧪 Generated {len(mutants)} mutants")

    mut_fasta = f"{Path(cleaned_fasta).stem}_mutants.fasta"
    write_fasta([(f"{sid}_{mid}", s) for sid, mid, s in mutants],
                mut_fasta)

    X = extract_features(mut_fasta)
    probs, preds = predict(model, X, threshold)

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
    print("# Prediction of Allergenic and Non-Allergenic Peptides                       #")
    print("# Raghava Group | IIIT-Delhi                                                 #")
    print("##############################################################################")

    args = parse_args()

    open(LOG_FILE, "w").close()

    # 🔥 AUTO SETUP
    setup_dependencies()

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


if __name__ == "__main__":
    main()
