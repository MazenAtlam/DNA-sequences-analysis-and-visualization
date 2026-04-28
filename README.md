# DNA Sequence Analysis and Visualization

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![Biopython](https://img.shields.io/badge/Biopython-Enabled-brightgreen)
![NumPy](https://img.shields.io/badge/NumPy-Data-blueviolet)
![Matplotlib](https://img.shields.io/badge/Matplotlib-Visuals-orange)
![License](https://img.shields.io/badge/License-MIT-green.svg)

A comprehensive Python toolkit for bioinformatics sequence analysis. This project provides programmatic modules for Open Reading Frame (ORF) detection, structural sequence alignment via 2D Dot Plots, and Codon Usage profiling.

Developed as part of the Bioinformatics curriculum at the Systems and Biomedical Engineering (SBME) department, Cairo University.

---

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Project Structure](#project-structure)
- [Quick Start](#quick-start)
- [Usage Guide](#usage-guide)
  - [1. ORF Detection](#1-orf-detection)
  - [2. Dot Plot Alignment](#2-dot-plot-alignment)
  - [3. Codon Usage Table](#3-codon-usage-table)
- [License](#license)

---

## Features

- **Open Reading Frame (ORF) Detection:** Safely scans DNA sequences across all three forward reading frames to identify valid ORFs (Start: `ATG`, Stop: `TAA`, `TAG`, `TGA`). Includes a customizable minimum nucleotide length filter.
- **Sequence Alignment Visualization:** Generates 2D dot plots to visualize exact character matches, inverted repeats, and indels between two DNA sequences using Matplotlib.
- **Codon Usage Profiling:** Calculates absolute codon frequencies for single DNA sequences and aggregates global codon usage distributions from multi-sequence `.fasta` files.

---

## Prerequisites

To run the scripts in this repository, you will need **Python 3.8 or higher** and the following external libraries:

- `biopython`: For secure biological sequence parsing and file I/O operations.
- `numpy`: For high-performance 2D matrix initialization and coordinate mapping.
- `matplotlib`: For rendering programmatic visualizations and charts.

You can install all dependencies via pip:

```bash
pip install biopython numpy matplotlib
```

---

## Project Structure

```text
dna-sequences-analysis-and-visualization/
├── Open Reading Frames (ORF) Detection/
│   ├── orf_detection.py
│   └── modified_orf_detection.py
├── Sequence Alignment Visualization/
│   └── get_dot_plot.py
├── Codon Usage Table/
│   ├── get_codon_freq_table.py
│   └── get_combined_codon_freq_table.py
├── datasets/
│   ├── exFasta.fasta
│   ├── exFasta3Rec.fasta
│   ├── exons.txt
│   └── genomic_dna.txt
├── .gitignore
├── LICENSE
└── README.md
```

---

## Quick Start

1. **Clone the repository:**

   ```bash
   git clone [https://github.com/mazenatlam/DNA-sequences-analysis-and-visualization.git](https://github.com/mazenatlam/DNA-sequences-analysis-and-visualization.git)
   cd DNA-sequences-analysis-and-visualization
   ```

2. **Run a demonstration module:**
   Each script contains a stress-testing block at the bottom of the file. You can execute them directly to see the output.

   ```bash
   # Example: Run the combined FASTA codon frequency generator
   python "Codon Usage Table/get_combined_codon_freq_table.py"
   ```

---

## Usage Guide

### 1. ORF Detection

Extract all Open Reading Frames from a sequence, or filter them by length.

```python
from Bio.Seq import Seq
import sys
sys.path.append('./Open Reading Frames (ORF) Detection')
from modified_orf_detection import get_orfs_with_min_length

my_seq = Seq("ATGCGTAAATGATAGTATAAAATGCCCTGA")
# Extract ORFs that are at least 10 nucleotides long
valid_orfs = get_orfs_with_min_length(my_seq, min_length=10)

for orf in valid_orfs:
    print(f"Frame {orf['frame']}: {orf['sequence']} (Length: {orf['length']})")
```

### 2. Dot Plot Alignment

Compare two sequences visually to find structural similarities.

```python
from Bio.Seq import Seq
import sys
sys.path.append('./Sequence Alignment Visualization')
from get_dot_plot import generate_dot_plot

seq_a = Seq("ATGCGTAAATGATAGTATAAAA")
seq_b = Seq("ATGATAGTAT")

# This will open a Matplotlib window rendering the 2D grid
generate_dot_plot(seq_a, seq_b)
```

### 3. Codon Usage Table

Generate a bar chart detailing the frequency of each codon from a standard FASTA dataset.

```python
import sys
sys.path.append('./Codon Usage Table')
from get_combined_codon_freq_table import render_combined_codon_table

# Parse the dataset and render the global frequency chart
render_combined_codon_table('./datasets/exFasta3Rec.fasta')
```

---

## License

This project is licensed under the terms of the MIT license.
