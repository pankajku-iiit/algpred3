# AlgPred3

AlgPred3 is a computational framework for the prediction, analysis, and design of **allergenic peptides** using machine learning, sequence-derived features, and motif-based approaches.

## Introduction

AlgPred3 is specifically designed for peptide-level allergen prediction using curated allergenic and non-allergenic peptide datasets. The final prediction model is based on an **Extra Trees (ET)** classifier trained using an optimized feature subset selected through **SVC-L1 feature selection**.

### Key Features
- Prediction of allergenic and non-allergenic peptides
- Protein scanning for allergenic regions
- Peptide design and mutant generation
- Motif-based allergen detection
- SHAP-based model interpretability
- Web server, standalone package, GitHub repository, and PyPI package

## Web Server
https://webs.iiitd.edu.in/raghava/algpred3/

## Installation

```bash
pip install algpred3
```

## GitHub Installation

```bash
git clone https://github.com/raghavagps/algpred3.git
cd algpred3
python algpred3.py -h
```

## Usage

### Prediction

```bash
python algpred3.py -i example.fasta -o prediction.csv -j pred
```

### Protein Scan

```bash
python algpred3.py -i protein.fasta -o scan.csv -j scan
```

### Design

```bash
python algpred3.py -i peptide.fasta -o design.csv -j des
```

## Workflow

1. Sequence validation
2. Feature extraction
3. Feature selection (SVC-L1)
4. Extra Trees prediction
5. SHAP interpretation
6. Output generation

## Citation

Kumar N, Mehta NK, Kumar P, Bajiya N, Rathore AS, Raghava GPS.

**AlgPred3: Prediction of Allergenic Peptides and Identification of Allergen-Specific Sequence Motifs.**

(Journal information to be added after publication)

## Disclaimer

AlgPred3 is intended for research purposes only. Predictions should be considered computational assessments and hypothesis-generating results. Experimental validation is recommended.
