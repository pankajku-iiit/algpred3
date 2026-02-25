AlgPred 3.0: Allergenicity Prediction & Peptide Design Suite

AlgPred 3.0 is a comprehensive toolkit for predicting allergenic and non-allergenic peptides using machine learning models and sequence-derived descriptors. It integrates automated feature extraction, prediction, protein scanning, and mutation design workflows in a unified framework.

It supports:

Standalone allergenicity prediction
Protein region scanning (sliding window analysis)
Exhaustive single mutation design analysis

The toolkit is optimized for Linux/macOS environments and supports reproducible command-line execution.

---

## Overview

AlgPred 3.0 predicts allergenic potential of peptide or protein sequences using:

Machine learning classifier trained on curated allergen datasets
Sequence descriptors generated via the **pfeature_comp** feature engine

It provides three workflows:

Prediction Mode — Predict allergenicity of peptides/proteins
Scan Mode — Identify allergenic regions in proteins
Design Mode — Generate and evaluate single mutants

---

## Core Functionalities

Allergenicity prediction from sequences
Automatic FASTA cleaning and validation
Sequence descriptor calculation via pfeature
Protein allergenicity scanning
Single mutation design analysis
Automated model download if absent

---

## Installation

### Option 1 — Conda Environment (Recommended)

```bash
conda create -n algpred3 python=3.10
conda activate algpred3
conda install pandas joblib scikit-learn
```

Clone repository:

```bash
git clone https://github.com/raghavagps/algpred3.git
cd algpred3
```

Download feature extraction binary:

```bash
wget <pfeature_release_link>
chmod +x pfeature_comp
```

---

### Option 2 — Manual Installation

Install dependencies:

```bash
pip install pandas joblib scikit-learn
```

Ensure:

* `pfeature_comp` executable is present
* `columns.csv` file exists
* Internet access for automatic model download

---

## Usage Overview

General syntax:

```bash
python algpred3.py -i input.fasta -j <mode> -o output.csv
```

Modes available:

* `pred` → Prediction mode
* `scan` → Protein scanning
* `des` → Mutation design

---

## Prediction Mode

Predict allergenicity of sequences.

```bash
python algpred3.py -i peptides.fasta -j pred -o output.csv
```

### Output Columns

Sequence_ID — Sequence identifier
Probability — Prediction score
Status — Allergen / Non-Allergen

---

## Scan Mode

Sliding window allergenicity analysis across proteins.

```bash
python algpred3.py -i protein.fasta -j scan -l 15 -s 1 -o scan.csv
```

Arguments:

-l → Window length (required)
-s → Step size (default: 1)

### Output Columns

ParentSeq — Protein ID
Start / End — Window positions
Peptide — Extracted fragment
Score — Prediction probability
Prediction — Allergen / Non-Allergen

---

## Design Mode

Generates all possible single amino acid mutants and evaluates allergenicity.

```bash
python algpred3.py -i peptides.fasta -j des -o design.csv
```

### Output Columns

SeqID — Original sequence ID
MutantID — Mutation annotation
Sequence — Mutant sequence
Score — Prediction probability
Prediction — Allergen / Non-Allergen

---

## Sequence Validation

Only standard amino acids allowed:

ACDEFGHIKLMNPQRSTVWY

Invalid sequences are:

Automatically removed
Logged in `stand_error.log`

Clean FASTA files are generated before analysis.

---

## Machine Learning Model

Model type:
Classification model (serialized Joblib format)

Feature source:
pfeature sequence descriptors

Model file:

`algpred3_model.sav`
Downloaded automatically if missing.

---

## Repository Contents

algpred3.py — Main unified script
algpred3_model.sav — Pretrained model (auto-download)
pfeature_comp — Descriptor generator binary
columns.csv — Selected feature columns
README.md — Documentation

---

## Citation

If you use AlgPred 3.0 in research, please cite the corresponding publication from:

Raghava GPS Group, IIIT-Delhi.

---

## Support

GitHub:
https://github.com/raghavagps/algpred3

Email:
raghava@iiitd.ac.in

---

## License

Academic and research use recommended.
Refer to the repository license file for details.
