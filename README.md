# 🧬 AlgPred 3.0

**Prediction of Allergenic and Non-Allergenic Peptides**
Developed by Raghava Group, IIIT-Delhi

---

## 📌 Overview

AlgPred 3.0 is a command-line tool for predicting whether peptide sequences are **allergenic or non-allergenic** using a trained machine learning model.

It supports three major modes:

* 🔍 **Prediction (`pred`)** – Classify full sequences
* 🧬 **Scan (`scan`)** – Sliding window analysis of proteins
* 🧪 **Design (`des`)** – Generate and evaluate mutants

---

## ⚙️ Features

* Automatic dependency download (model + binaries)
* FASTA & plain text input support
* Built-in sequence cleaning and validation
* Multiple feature extraction techniques:

  * AAC, DPC, BTC, DDR
  * PAAC, APAAC, CeTD
* Adjustable prediction threshold
* CSV output for easy analysis

---

## 📂 Project Structure

```
.
├── algpred3.py
├── algpred3_model.sav      # Auto-downloaded
├── pfeature_comp           # Auto-downloaded
├── columns.csv             # Auto-downloaded
├── stand_error.log         # Invalid sequences log
└── output files (.csv)
```

---

## 🚀 Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-username/algpred3.git
cd algpred3
```

### 2. Install dependencies

```bash
pip install pandas joblib
```

---

## ▶️ Usage

```bash
python algpred3.py -i INPUT -j JOB [options]
```

### 🔹 Required Arguments

| Argument | Description                      |
| -------- | -------------------------------- |
| `-i`     | Input file (FASTA or plain text) |
| `-j`     | Job type: `pred`, `scan`, `des`  |

---

### 🔹 Optional Arguments

| Argument | Default               | Description               |
| -------- | --------------------- | ------------------------- |
| `-o`     | final_predictions.csv | Output file               |
| `-t`     | 0.5                   | Prediction threshold      |
| `-l`     | —                     | Window length (scan mode) |
| `-s`     | 1                     | Step size (scan mode)     |

---

## 🧪 Modes Explained

### 🔍 1. Prediction Mode (`pred`)

Predict allergenicity of full sequences:

```bash
python algpred3.py -i input.fasta -j pred
```

Output:

* Sequence ID
* Probability score
* Allergen / Non-Allergen

---

### 🧬 2. Scan Mode (`scan`)

Perform sliding window analysis:

```bash
python algpred3.py -i input.fasta -j scan -l 20 -s 5
```

Output:

* Parent sequence
* Start & end positions
* Peptide fragment
* Prediction

---

### 🧪 3. Design Mode (`des`)

Generate all possible single-point mutants:

```bash
python algpred3.py -i input.fasta -j des
```

Output:

* Mutant ID
* Modified sequence
* Prediction score

---

## 🧹 Input Handling

* Accepts:

  * FASTA format
  * Plain text sequences
* Automatically:

  * Removes invalid amino acids
  * Logs removed sequences in `stand_error.log`

Valid amino acids:

```
ACDEFGHIKLMNPQRSTVWY
```

---

## 📊 Output

All results are saved in **CSV format**, including:

* Prediction scores
* Classification labels
* Sequence metadata (depending on mode)

---

## 📦 Auto-Downloaded Dependencies

The script automatically downloads:

* Trained ML model
* Feature extraction binary (`pfeature`)
* Required feature column list

---

## ⚠️ Notes

* Ensure internet connection for first run (auto-downloads)
* `scan` mode requires `-l` (window length)
* Large sequences in `des` mode may generate **many mutants**

---

## 🙏 Acknowledgment

Developed by **Raghava Group**,
Indraprastha Institute of Information Technology, Delhi (IIIT-Delhi)

---

## 💡 Example

```bash
python algpred3.py -i sample.fasta -j pred -t 0.6 -o results.csv
```
