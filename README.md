# AlgPred3

A computational framework for predicting and designing **allergenic and non-allergenic peptides** using machine learning and composition-based features.

---

## 📌 Introduction

**AlgPred3** is developed to identify allergenic potential of peptides and proteins based on their primary sequence. It integrates feature-based machine learning approaches with curated allergen datasets to provide robust predictions.
🔗 Visit the web server for more information: [Algpred3](https://webs.iiitd.edu.in/raghava/algpred3/)
The tool supports:

* **Prediction of allergenic peptides**
* **Protein scanning for allergenic regions**
* **Design of peptide mutants**

It uses compositional and physicochemical descriptors such as:

* AAC (Amino Acid Composition)
* DPC (Dipeptide Composition)
* PAAC, APAAC
* CeTD, DDR, BTC

---

## 📚 Reference

**AlgPred3**
Raghava Group, IIIT-Delhi

---

## 🧪 Quick Start for Reproducibility

```bash
# 1. Clone the repository
git clone https://github.com/raghavagps/algpred3.git
cd algpred3
# 2. Set up the environment (conda recommended)
conda env create -f environment.yml
conda activate algpred3

# 3. Run help
python algpred3.py -h

# 4. Run prediction
python algpred3.py -i example.fasta -o output.csv -j pred
```

---

## 🛠️ Installation Options

### 🧰 Pip Installation

```bash
pip install algpred3
```

Check options:

```bash
algpred3 -h
```

---

### 🔹 Standalone Installation

AlgPred3 is written in **Python 3** and requires:

#### ✅ Required Libraries

```bash
python=3.8+
```

#### 📦 Install Dependencies

```bash
pip install scikit-learn
pip install pandas
pip install joblib
```

---

### 🔹 Automatic Dependency Setup

AlgPred3 automatically downloads required files during runtime:

* Model file (`algpred3_model.sav`)
* Feature column file (`columns.csv`)
* pfeature binary

This behavior is implemented directly in the script. 

---

## ⚠️ Important Notes

* First run requires **internet connection** for auto-download.
* Input sequences must contain only valid amino acids:
  **ACDEFGHIKLMNPQRSTVWY**
* Invalid sequences are removed and logged in:

```bash
stand_error.log
```

---

## 🔬 Classification

AlgPred3 classifies peptides into:

* **Allergen**
* **Non-Allergen**

based on machine learning predictions using extracted sequence features.

---

## 🚀 Usage

### 🔹 Minimum Usage

```bash
python algpred3.py -h
```

---

### 🔹 Full Usage

```bash
usage: algpred3.py [-h]
                  -i INPUT
                  [-o OUTPUT]
                  -j {pred,scan,des}
                  [-t THRESHOLD]
                  [-l LENGTH]
                  [-s STEP]
```

---

### 📌 Arguments

| Argument    | Description                                        |
| ----------- | -------------------------------------------------- |
| `-i INPUT`  | Input FASTA or sequence file                       |
| `-o OUTPUT` | Output CSV file (default: `final_predictions.csv`) |
| `-j`        | Job type: `pred`, `scan`, `des`                    |
| `-t`        | Threshold (default: 0.5)                           |
| `-l`        | Window length (scan mode only)                     |
| `-s`        | Step size (default: 1)                             |

---

## 📂 Input & Output Files

### ✅ Input Formats

1. **FASTA format**
2. **Plain text format** (one sequence per line)

---

### ✅ Output File

* Results are saved in **CSV format**
* Includes:

  * Sequence ID
  * Prediction score
  * Allergen / Non-Allergen status

---

## 🔍 Jobs & Features

### 🔹 Job Types

| Job         | Description                                |
| ----------- | ------------------------------------------ |
| 🔹 **pred** | Predict allergenic potential of sequences  |
| 🔹 **scan** | Identify allergenic regions in proteins    |
| 🔹 **des**  | Generate and evaluate all possible mutants |

---

### 🔹 Scan Mode

* Uses **sliding window approach**
* Generates peptide fragments
* Predicts allergenic regions

---

### 🔹 Design Mode

* Generates **all single amino acid mutants**
* Predicts allergenicity of each mutant
* Useful for:

  * Peptide optimization
  * Allergenicity reduction

---

## ⚙️ Workflow

1. Input sequence validation
2. Cleaning invalid sequences
3. Feature extraction (pfeature)
4. Feature selection
5. Model prediction
6. Output generation

---

## 📑 Package Contents

| File                   | Description               |
| ---------------------- | ------------------------- |
| **algpred3.py**        | Main prediction script    |
| **columns.csv**        | Feature selection file    |
| **algpred3_model.sav** | Trained ML model          |
| **pfeature_comp**      | Feature extraction binary |
| **example.fasta**      | Sample input              |

---

## 📦 PIP Installation (Again for Reference)

```bash
pip install algpred3
```

---

🚀 **Start predicting allergenic peptides with AlgPred3 today!**
🔗 Visit: [Algpred3](https://webs.iiitd.edu.in/raghava/algpred3/)
---
