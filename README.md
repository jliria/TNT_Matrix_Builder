# TNT Matrix Builder: A Python Graphical Interface for Total Evidence Phylogenetics

<div align="left">
  <img src="https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEhEhENEvuyWLIKZ_QziPGKwxQOYWxh7bT8u2Y0v-MMvjORJ1wQHt_-8H6x7z7IGbjyuRCL5BFLjmV2DyJmmw6EFwjI08TH1jFl0G4WyMGmMpYq6qOBhVj5QSyDHNgqtHhnx7182Mg-T2lNjIXkTagXYdlHo4SE7yO2isrnGP3UM0fSt2MoPbSA3DN8LHk6k/w286-h236/Untitled.jpg" width="250" alt="TNT Matrix Builder Logo" style="border-radius: 8px; margin-bottom: 20px;" />
 
  [![Python Version](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
</div>

---

## Introduction

`TNT Matrix Builder` is a graphical application developed in Python specifically designed to facilitate the teaching and learning of the TNT software in Phylogenetic Systematics and Cladistics courses. The program addresses one of the greatest challenges for students and researchers: the integration of "total evidence" (molecular data, discrete morphological data, and continuous characters) into a single analysis file. 

While packages like Mesquite (Maddison & Maddison, 2025) are excellent for the visual editing of multiple matrices, they lack a combined export function that automates the syntax required by TNT[cite: 1]. This forces the user to resort to manual text editing, a process involving the precise configuration of commands like `nstates`, `xread`, and block tags (`&[dna]`, `&[num]`, `&[cont]`), which are frequently sources of errors during practical exercises.

This tool complements and expands the workflow of previous applications such as `py_tm2tnt` and `py_tps2tnt` (Liria & Soto-Vivas, 2025), allowing students or researchers who are new to TNT to perform phylogenetic analyses incorporating an integrative approach. By automating the concatenation of multiple data blocks, `TNT Matrix Builder` eliminates potential errors in the manipulation of the `.tnt` file, allowing the user to focus strictly on the phylogenetic analysis as such. 

The program not only automatically synchronizes taxon names using union and intersection protocols, but also manages the complexity of 2D geometric morphometric data (TPS files), averaging specimens per taxon and integrating measurement scales according to the guidelines of Catalano & Goloboff (2018). This tool is therefore presented as an essential technological bridge that simplifies the transition from heterogeneous data collection to its processing in the TNT phylogenetic search engine (Goloboff & Morales, 2023).

---

## User Interface Overview

<div align="center">
  <!-- Replace with actual path to image if hosted locally -->
  <img src="images/figure1.png" width="700" alt="TNT Matrix Builder Interface Overview" style="border-radius: 6px; border: 1px solid #ddd;" />
  <p><em>Figure 1. General overview of the TNT Matrix Builder interface showing the "Blocks" area and the lateral configuration panel[cite: 1].</em></p>
</div>

---

## 1. Requirements

Before using `TNT Matrix Builder`, ensure you have installed[cite: 1]:
* **Python 3.x**[cite: 1]
* **Standard libraries:** `tkinter`, `re`, `math`, `os`, `dataclasses`[cite: 1].
* **Optional library:** `tkinterdnd2` (to enable the drag-and-drop file feature)[cite: 1]. You can install it using the command: `pip install tkinterdnd2`[cite: 1].

*Note: There is an executable version available that can be downloaded from the following link: [Insert Link Here]*[cite: 1].

---

## 2. Supported Formats

The program automatically detects the nature of the data based on its content[cite: 1]:
* **Molecular (DNA):** FASTA files (`.fasta`, `.fa`, `.fas`)[cite: 1].
* **Tables (`.txt`, `.tsv`, `.csv`):**
  * **Discrete/Numeric:** If the file contains integers or polymorphisms[cite: 1].
  * **Continuous:** Automatically detected by the presence of decimal values[cite: 1].
* **Morphometrics (TPS):** `.tps` files with 2D landmarks[cite: 1]. The program automatically averages the specimens per species[cite: 1].

<div align="center">
  <img src="images/figure2.png" width="700" alt="Supported Formats Example" style="border-radius: 6px; border: 1px solid #ddd;" />
  <p><em>Figure 2. Example of the various supported data types: FASTA sequences, continuous/discrete data tables, and morphometric TPS files[cite: 1].</em></p>
</div>

---

## 3. Step-by-Step Guide

### Step 1: Loading files
Add your matrices using the **"Add files..."** button or by dragging the files directly onto the top panel of the interface. The blocks will appear listed showing their automatically detected data type.

### Step 2: Organization and priority rule
This option allows organizing the order of the blocks in TNT; for example, if it is desired that molecular characters (COI, 16s, etc.) be at the end, and others (continuous, discrete) at the beginning. If the project includes continuous characters or Landmarks, **these must always be placed in the first position of the list (Block 1)**. 

Utilize the **"Move up"** and **"Move down"** buttons to adjust the sorting order. You can use **"Rename block..."** to assign descriptive names for each block, which will then appear in the `cname` section of the final file.

<div align="center">
  <img src="images/figure3.png" width="700" alt="Block Priority Sorting" style="border-radius: 6px; border: 1px solid #ddd;" />
  <p><em>Figure 3. Detail of the block list indicating the use of movement buttons to position continuous data at the beginning[cite: 1].</em></p>
</div>

### Step 3: Taxa concatenation (Taxa Merge)
* **Union:** Keeps all taxon names and fills missing blocks with empty/missing data (`?`)[cite: 1].
* **Intersection:** Only retains the taxa that are present across all loaded files[cite: 1].
* **Normalize names:** Replaces spaces with underscores (`_`) for terminal names where Genus and species are assigned, which is a strict requirement in TNT and a common source of formatting errors[cite: 1].

<div align="center">
  <img src="images/figure4.png" width="700" alt="Taxa Merge Options" style="border-radius: 6px; border: 1px solid #ddd;" />
  <p><em>Figure 4. Option panel for "Taxa merge" with the Union/Intersection functions and species name normalization[cite: 1].</em></p>
</div>

### Step 4: Exporting the matrix to TNT
Select the range of states in **"nstates"** (the program will automatically choose `32` if continuous data is detected)[cite: 1]. If the concatenated matrices only contain molecular data, it is highly recommended to choose `nstates 8` to optimize engine memory usage[cite: 1]. Finally, click on **"Export .tnt..."** to obtain your ready-to-run concatenated file[cite: 1].

---

## 4. Concatenation report

Before exporting, it is highly recommended to check the **"Show taxa report..."** button[cite: 1]. This function allows you to validate the structural integrity of the matrix, showing which species are absent in each individual block and confirming the total number of characters that will make up the final analysis[cite: 1].

<div align="center">
  <img src="images/figure5.png" width="700" alt="Taxa Report Window" style="border-radius: 6px; border: 1px solid #ddd;" />
  <p><em>Figure 5. Taxa report window ("Taxa report") showing the list of absent species for each loaded source file[cite: 1].</em></p>
</div>

---

## 5. References

* Catalano, S., Goloboff, P.A. (2018). *A guide for the analysis of continuous and landmark characters in TNT*.
* Goloboff, P.A., Morales, M.E. (2023). *TNT version 1.6, with a graphical interface for MacOS and Linux, including new routines in parallel.* Cladistics, 39: 144-153.
* Liria, J., Soto-Vivas, A. (2025). *py_tps2tnt y py_tm2tnt: dos programas en Python para procesamiento de datos morfométricos en análisis cladísticos con TNT.* Revista Peruana de Biología, 32(2), e30018.
* Maddison, W. P., D.R. Maddison. (2025). *Mesquite: a modular system for evolutionary analysis.* Version 4.02. https://www.mesquiteproject.org.

---

## License

This project is licensed under the MIT License - see the local repository files for details.
