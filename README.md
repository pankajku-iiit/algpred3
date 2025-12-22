# AlgPred 3.0

## Overview
AlgPred 3.0 is a standalone command-line tool for the **prediction, scanning, and design of allergenic and non-allergenic proteins/peptides** from their primary amino acid sequences. The tool is developed by the **Raghava Group**, Department of Computational Biology, **Indraprastha Institute of Information Technology, Delhi (IIIT-Delhi)**.

AlgPred 3.0 supports three major functionalities:

- **Prediction (pred):** Classify full-length protein/peptide sequences as allergen or non-allergen.
- **Protein Scan (scan):** Sliding-window scanning of proteins to identify allergenic regions.
- **Design (des):** Systematic single-point mutation analysis to design reduced-allergenicity variants.

---

## Features

- Accepts FASTA or plain-text sequence input
- Automatic sequence cleaning and validation
- Dipeptide Composition (DPC) based feature extraction
- Machine-learning-based probability prediction
- Sliding window allergen mapping (Protein Scan)
- Exhaustive single-point mutant generation (Design mode)
- CSV outputs suitable for web-server integration
- Compatible with PHP/HTML visualization pipelines

---

## Requirements

- Python 3.7 or higher
- Required Python libraries:
  - pandas
  - joblib

Install dependencies using:

```bash
pip install pandas joblib
```

---

## Model File

MOdel file will be atomatically downloaded from the website of AlgPred 3.0. Model can be manually downloaded from 
https://webs.iiitd.edu.in/raghava/algpred3/algpred3_model.sav
```
algpred3_model.sav
```

The program will terminate if the model file is missing.

---

## Input Format

### FASTA format

```text
>seq1
MAVPQNRVTRSRRNMRRAHDALVAANPASCPNCGELKRPHHVCGACGHYDDREVVAQAAEVDLDDDAA
>seq2
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQLR
```

### Plain text format

```text
MAVPQNRVTRSRRNMRRAHDALVAANPASCPNCGELKRP
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQLR
```

Invalid characters (non-standard amino acids) are automatically removed and logged.

---

## Usage

```bash
python algpred3.py -i input.fasta -j JOBTYPE [options]
```

### Common arguments

| Argument | Description |
|--------|-------------|
| `-i, --input` | Input FASTA or plain text file |
| `-o, --output` | Output CSV file (default: final_predictions.csv) |
| `-j, --job` | Job type: pred, map, or des |
| `-t, --threshold` | Prediction threshold (default: 0.5) |
| `-wd, --working_directory` | Working directory |

---

## Job Modes

### 1. Prediction Mode (`pred`)

Predict allergenicity of full-length sequences.

```bash
python3 algpred3.py -i input.fasta -j pred -o pred_output.csv
```

**Output columns:**
- Sequence_ID
- Probability
- Status (Allergen / Non-Allergen)

---

### 2. Protein Scan Mode (`scan)

Performs sliding-window allergen scanning across protein sequences.

```bash
python3 algpred3.py -i input.fasta -j scan -l 10 -s 1 -o scan_output.csv
```

**Required arguments:**
- `-l, --length` : Sliding window length

**Optional:**
- `-s, --step` : Step size (default: 1)

**Output columns:**
- ParentSeq
- Start
- End
- Peptide
- Score
- Prediction

This output is designed for downstream residue-level mapping and visualization.

---

### 3. Design Mode (`des`)

Generates all possible single-point mutants and predicts allergenicity.

```bash
python3 algpred3.py -i input.fasta -j des -o design_output.csv
```

**Output columns:**
- SeqID
- MutantID
- Sequence
- Score
- Prediction

Design mode can generate very large outputs for long sequences.

---

## Logs

- Invalid sequences and removed residues are logged to:

```
stand_error.log
```

---

## Output Files

- Cleaned FASTA: `*_clean.fasta`
- Feature files: `*_DPC.csv`
- Scan windows FASTA: `scan_windows.fasta`
- Final prediction CSVs as specified by `-o`

---

## Citation

If you use AlgPred 3.0 in your research, please cite:

> **AlgPred 3.0**  
> Raghava Group, Department of Computational Biology  
> Indraprastha Institute of Information Technology, Delhi (IIIT-Delhi)

(Please update this section with journal reference and DOI when published.)

---

## License

This software is intended for **academic and research use**. Redistribution or commercial use requires permission from the authors.

---

## Contact

**Raghava Group**  
Department of Computational Biology  
IIIT-Delhi, India

Website: https://webs.iiitd.edu.in/raghava/

---

## Acknowledgements

We acknowledge the developers and users of AlgPred series for continuous feedback and improvements.

---

Thank you for using **AlgPred 3.0**.

