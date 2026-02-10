# TNT_Matrix_Builder
TNT Matrix Builder is a specialized graphical application designed to streamline the construction of "Total Evidence" matrices for the phylogenetic software TNT (Tree Analysis using New Technologies).

Building complex matrices often requires manual editing of .tnt files to include specific headers (nstates, xread) and block tags (&[dna], &[num], &[cont]). This tool automates the entire process, allowing researchers and students to focus on biological hypotheses rather than syntax debugging.

Key Features:
Multi-Format Support: Automatically detects and parses FASTA (DNA), TPS (2D Landmarks), and tabular data (CSV/TSV/TXT for discrete or continuous characters).

Smart Taxa Merging: Handles disparate datasets using Union or Intersection modes, automatically filling gaps with missing data (?).

Morphometric Integration: Built-in support for 2D landmarks, including specimen averaging and scale management based on Catalano & Goloboff (2018) guidelines.

Validation Tools: Provides a "Taxa Report" to identify naming inconsistencies or missing species across blocks before exporting.

Automatic TNT Syntax: Generates perfectly formatted .tnt files, including the cname section for easy character identification.

Technical Stack:
Language: Python 3.x

GUI: Tkinter (with optional tkinterdnd2 support for drag-and-drop).

Dependencies: Standard library (math, re, os, dataclasses).
