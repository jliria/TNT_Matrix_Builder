"""
TNT Matrix Builder (GUI) - COMPLETE v1 (with TPS landmarks 2D) - Liria & Soto-Vivas (2026)

Features:
- Drag & drop multiple files (optional) or Add Files...
- Detects data type: DNA (FASTA), continuous table, discrete/numeric table, TPS landmarks 2D
- Validates & merges taxa (Union or Intersection)
- Orders blocks, names blocks, and exports a .tnt file compatible with TNT
- Designed to be improved iteratively during testing

Supported:
- FASTA: .fasta .fa .fas  -> &[dna]
- Tables: .txt .tsv .csv   -> auto-detect &[cont] vs &[num]
- TPS:   .tps              -> &[landmark 2D]  (averages multiple specimens per species)

Optional drag&drop:
    pip install tkinterdnd2
"""

from __future__ import annotations

import os
import re
import math
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Any

HAS_DND = False
try:
    from tkinterdnd2 import DND_FILES, TkinterDnD 
    HAS_DND = True
except Exception:
    HAS_DND = False


# -------------------------
# Utilities
# -------------------------
def sanitize_taxon(name: str) -> str:
    """Normalize taxon names for TNT: spaces -> '_' and remove weird chars."""
    s = name.strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^\w\-\.\:]", "_", s)
    return s


def split_dnd_files(raw: str) -> List[str]:
    """
    Split drag&drop paths.
    Handles: {C:\My File.txt} C:\Other.txt
    """
    out, buf, in_brace = [], "", False
    for ch in raw:
        if ch == "{":
            in_brace = True
            buf = ""
        elif ch == "}":
            in_brace = False
            if buf.strip():
                out.append(buf.strip())
            buf = ""
        elif ch == " " and not in_brace:
            if buf.strip():
                out.append(buf.strip())
                buf = ""
        else:
            buf += ch
    if buf.strip():
        out.append(buf.strip())
    return out


def sniff_delimiter(first_line: str) -> str:
    """Rudimentary delimiter sniff for CSV/TSV."""
    if "\t" in first_line:
        return "\t"
    if "," in first_line:
        return ","    
    return " "

def is_missing(x: str) -> bool:
    x = x.strip()
    return x == "" or x.upper() in {"NA", "N/A"} or x in {"?", "-"}

def is_int_like(s: str) -> bool:
    s = s.strip()
    if s.startswith("{") and s.endswith("}"):
        # TNT polymorphism like {01}
        inner = s[1:-1].strip()
        return bool(re.fullmatch(r"[0-9]+", inner))
    return bool(re.fullmatch(r"[0-9]+", s))

def parse_float_safe(s: str) -> Optional[float]:
    try:
        return float(s)
    except Exception:
        return None

# -------------------------
# Data structures
# -------------------------
@dataclass
class DataBlock:
    source_path: str
    name: str
    dtype: str  # 'dna' | 'cont' | 'num' | 'landmark2d'
    taxa_original: List[str]
    taxa: List[str]  
    nchar: int
    rows_by_taxon: Dict[str, Any]  # taxon -> row
    missing_symbol: str = "?"
    warnings: List[str] = field(default_factory=list)
    meta: Dict[str, Any] = field(default_factory=dict)

    def short_type(self) -> str:
        return self.dtype.upper()

# -------------------------
# Parsers / Detection
# -------------------------
class Detector:
    @staticmethod
    def detect(path: str) -> str:
        ext = os.path.splitext(path.lower())[1]
        if ext in [".fasta", ".fa", ".fas"]:
            return "fasta"
        if ext in [".csv", ".tsv", ".txt"]:
            return "table"
        if ext == ".tps":
            return "tps"
        return "unknown"

class Parsers:
    @staticmethod
    def parse(path: str) -> DataBlock:
        kind = Detector.detect(path)
        base_name = os.path.splitext(os.path.basename(path))[0]

        if kind == "fasta":
            return Parsers.parse_fasta(path, base_name)
        if kind == "table":
            return Parsers.parse_table(path, base_name)
        if kind == "tps":
            return Parsers.parse_tps_2d(path, base_name)
        raise ValueError(f"Unsupported/unknown file type: {path}")

    @staticmethod
    def parse_fasta(path: str, name: str) -> DataBlock:
        taxa_original: List[str] = []
        seqs: Dict[str, str] = {}

        current = None
        chunks: List[str] = []

        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current is not None:
                        seqs[current] = "".join(chunks).upper()
                    current = line[1:].strip()
                    taxa_original.append(current)
                    chunks = []
                else:
                    chunks.append(re.sub(r"\s+", "", line))

        if current is not None:
            seqs[current] = "".join(chunks).upper()

        if not seqs:
            raise ValueError("FASTA appears empty or invalid.")

        # Determine alignment length
        lengths = {len(s) for s in seqs.values()}
        if len(lengths) != 1:
            raise ValueError(f"FASTA not aligned (multiple lengths found): {sorted(lengths)}")
        nchar = lengths.pop()

        taxa_sanitized = [sanitize_taxon(t) for t in taxa_original]
        if len(set(taxa_sanitized)) != len(taxa_sanitized):
            raise ValueError("FASTA taxon names collide after sanitization. Please rename taxa.")

        rows_by_taxon = {sanitize_taxon(t): seqs[t] for t in taxa_original}

        return DataBlock(
            source_path=path,
            name=name,
            dtype="dna",
            taxa_original=taxa_original,
            taxa=taxa_sanitized,
            nchar=nchar,
            rows_by_taxon=rows_by_taxon,
            missing_symbol="?",
            warnings=[]
        )

    @staticmethod
    def parse_table(path: str, name: str) -> DataBlock:
        """
        Reads a TSV/CSV/TXT matrix:
            taxon <delim> char1 <delim> char2 ...
        Auto-detect cont vs num:
            - cont if most non-missing cells parse as float AND at least one has decimal or non-int float
            - else num (discrete/numeric)
        """
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            lines = [ln.rstrip("\n") for ln in f if ln.strip()]

        if not lines:
            raise ValueError("Table file appears empty.")

        delim = sniff_delimiter(lines[0])
        taxa_original: List[str] = []
        data_rows: List[List[str]] = []

        for ln in lines:
            if delim == " ":
                parts = ln.split()
            else:
                parts = ln.split(delim)

            # Trim trailing empty columns (common in TSV)
            while parts and parts[-1].strip() == "":
                parts.pop()

            if len(parts) < 2:
                continue

            tax = parts[0].strip()
            vals = [p.strip() for p in parts[1:]]

            taxa_original.append(tax)
            data_rows.append(vals)

        if not data_rows:
            raise ValueError("No valid rows found in table.")

        nchar_set = {len(r) for r in data_rows}
        if len(nchar_set) != 1:
            raise ValueError(f"Inconsistent number of columns across rows: {sorted(nchar_set)}")
        nchar = nchar_set.pop()

        taxa_sanitized = [sanitize_taxon(t) for t in taxa_original]
        if len(set(taxa_sanitized)) != len(taxa_sanitized):
            raise ValueError("Table taxon names collide after sanitization. Please rename taxa.")

        # Detect cont vs num
        non_missing: List[str] = []
        for r in data_rows:
            for v in r:
                if not is_missing(v):
                    non_missing.append(v)

        dtype = "num"
        warnings: List[str] = []

        if non_missing:
            float_ok = 0
            has_decimal = False
            has_nonint_float = False
            for v in non_missing[:5000]:
                fv = parse_float_safe(v)
                if fv is not None:
                    float_ok += 1
                    if "." in v or "e" in v.lower():
                        has_decimal = True
                    if not is_int_like(v):
                        has_nonint_float = True
            frac_float = float_ok / max(1, len(non_missing))

            if frac_float >= 0.85 and (has_decimal or has_nonint_float):
                dtype = "cont"
            else:
                dtype = "num"

        rows_by_taxon: Dict[str, Any] = {}
        for tax, row in zip(taxa_sanitized, data_rows):
            rows_by_taxon[tax] = row

        return DataBlock(
            source_path=path,
            name=name,
            dtype=dtype,
            taxa_original=taxa_original,
            taxa=taxa_sanitized,
            nchar=nchar,
            rows_by_taxon=rows_by_taxon,
            missing_symbol="?",
            warnings=warnings,
        )

    @staticmethod
    def parse_tps_2d(path: str, name: str) -> DataBlock:
        """
        Parse TPS landmarks 2D.
        Expected (common) structure:
            LM=K
            x y   (K lines)
            ...
            ID=Taxon Name
        Repeats for multiple specimens.
        If a taxon appears multiple times, averages coordinates (landmark-wise).

        Also supports (best-effort):
            SCALE=val  (applies to the specimen where it appears)
        """
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            raw_lines = [ln.strip() for ln in f if ln.strip()]

        if not raw_lines:
            raise ValueError("TPS file appears empty.")

        specimens: List[Tuple[str, List[Tuple[float, float]]]] = []
        current_k: Optional[int] = None
        current_coords: List[Tuple[float, float]] = []
        current_id: Optional[str] = None
        current_scale: float = 1.0

        def flush_specimen():
            nonlocal current_k, current_coords, current_id, current_scale
            if current_k is None:
                return
            if len(current_coords) != current_k:
                raise ValueError(f"TPS specimen has LM={current_k} but read {len(current_coords)} coords.")
            if current_id is None:
                raise ValueError("TPS specimen missing ID=... line.")
            # apply scale
            if not math.isclose(current_scale, 1.0):
                current_coords = [(x * current_scale, y * current_scale) for x, y in current_coords]
            specimens.append((current_id, current_coords))
            # reset
            current_k = None
            current_coords = []
            current_id = None
            current_scale = 1.0

        for ln in raw_lines:
            up = ln.upper()

            if up.startswith("LM="):
                # flush previous specimen if partially accumulated
                if current_k is not None:
                    # If LM appears again, previous specimen should be complete; flush now.
                    flush_specimen()

                k_str = ln.split("=", 1)[1].strip()
                k = int(k_str)
                current_k = k
                current_coords = []
                current_id = None
                current_scale = 1.0
                continue

            if up.startswith("SCALE="):
                # applies to current specimen (or upcoming)
                s_str = ln.split("=", 1)[1].strip()
                try:
                    current_scale = float(s_str)
                except Exception:
                    raise ValueError(f"Invalid SCALE value: {ln}")
                continue

            if up.startswith("ID="):
                current_id = ln.split("=", 1)[1].strip()
                # once ID arrives, specimen is complete (in your file it comes after coords)
                flush_specimen()
                continue

            # coordinate line? (two numbers)
            if current_k is not None and len(current_coords) < current_k:
                parts = ln.split()
                if len(parts) >= 2:
                    x = parse_float_safe(parts[0])
                    y = parse_float_safe(parts[1])
                    if x is not None and y is not None:
                        current_coords.append((x, y))
                        continue

            # Ignore other TPS keys like IMAGE=, COMMENT=, etc.

        # flush if file ended mid-specimen (shouldn't, but handle)
        if current_k is not None:
            flush_specimen()

        if not specimens:
            raise ValueError("No specimens detected in TPS (need LM=... and ID=...).")

        # Group specimens by taxon and average if needed
        by_taxon: Dict[str, List[List[Tuple[float, float]]]] = {}
        for taxon, coords in specimens:
            by_taxon.setdefault(taxon, []).append(coords)

        taxa_original = list(by_taxon.keys())
        taxa_sanitized = [sanitize_taxon(t) for t in taxa_original]
        if len(set(taxa_sanitized)) != len(taxa_sanitized):
            raise ValueError("TPS taxon names collide after sanitization. Please rename taxa.")

        # Check LM count consistency
        lm_counts = {len(config[0]) for config in by_taxon.values()}
        if len(lm_counts) != 1:
            raise ValueError(f"TPS has inconsistent LM counts across taxa: {sorted(lm_counts)}")
        k = lm_counts.pop()

        rows_by_taxon: Dict[str, Any] = {}
        warnings: List[str] = []
        multi = [t for t, configs in by_taxon.items() if len(configs) > 1]
        if multi:
            warnings.append(f"{len(multi)} taxa had >1 specimen; used landmark-wise mean.")

        # Compute mean configuration per taxon
        for orig_tax, sanitized in zip(taxa_original, taxa_sanitized):
            configs = by_taxon[orig_tax]  # list of configurations, each length k
            if len(configs) == 1:
                mean_coords = configs[0]
            else:
                # landmark-wise mean
                mean_coords = []
                for i in range(k):
                    xs = [cfg[i][0] for cfg in configs]
                    ys = [cfg[i][1] for cfg in configs]
                    mean_coords.append((sum(xs) / len(xs), sum(ys) / len(ys)))
            rows_by_taxon[sanitized] = mean_coords

        return DataBlock(
            source_path=path,
            name=name,
            dtype="landmark2d",
            taxa_original=taxa_original,
            taxa=taxa_sanitized,
            nchar=1,  # landmarks treated as ONE character block in TNT, but with 2*k numbers per taxon
            rows_by_taxon=rows_by_taxon,
            missing_symbol="?",
            warnings=warnings,
            meta={"n_landmarks": k},
        )


# -------------------------
# Merge + Export
# -------------------------
class Project:
    def __init__(self):
        self.blocks: List[DataBlock] = []
        self.union_mode: str = "union"  # 'union' or 'intersection'
        self.normalize_names: bool = True
        self.decimals: int = 6
        self.nstates_mode: str = "auto"  # auto|8|16|32

        self.master_taxa: List[str] = []
        self.report: Dict[str, Dict[str, Any]] = {}

    def add_block(self, block: DataBlock):
        self.blocks.append(block)
        self.rebuild()

    def remove_block(self, idx: int):
        if 0 <= idx < len(self.blocks):
            self.blocks.pop(idx)
        self.rebuild()

    def move_block(self, idx: int, direction: int):
        j = idx + direction
        if 0 <= idx < len(self.blocks) and 0 <= j < len(self.blocks):
            self.blocks[idx], self.blocks[j] = self.blocks[j], self.blocks[idx]
        self.rebuild()

    def rebuild(self):
        self.build_master_taxa()
        self.validate()

    def build_master_taxa(self):
        if not self.blocks:
            self.master_taxa = []
            return

        sets = [set(b.taxa) for b in self.blocks]
        if self.union_mode == "intersection":
            master = set.intersection(*sets)
        else:
            master = set.union(*sets)

        self.master_taxa = sorted(master)

    def validate(self):
        rep: Dict[str, Dict[str, Any]] = {}
        master = set(self.master_taxa)

        for b in self.blocks:
            bt = set(b.taxa)
            rep[b.name] = {
                "dtype": b.dtype,
                "nchar": b.nchar,
                "ntax": len(b.taxa),
                "missing_taxa": sorted(master - bt),
                "extra_taxa": sorted(bt - master),
                "warnings": list(b.warnings),
                "source": b.source_path,
                "meta": dict(b.meta),
            }
        self.report = rep

    def total_nchar(self) -> int:
        return sum(b.nchar for b in self.blocks)

    def choose_nstates(self) -> int:
        if self.nstates_mode in {"8", "16", "32"}:
            return int(self.nstates_mode)
        # auto:
        has_cont = any(b.dtype == "cont" for b in self.blocks)
        if has_cont:
            return 32
        return 8

    def export_tnt(self, out_path: str, title: str = "Data saved from TNT"):
        if not self.blocks:
            raise ValueError("No blocks to export.")
        if not self.master_taxa:
            raise ValueError("Master taxa set is empty (check Union/Intersection).")

        nstates = self.choose_nstates()
        ntax = len(self.master_taxa)
        nchar_total = self.total_nchar()

        lines: List[str] = []
        lines.append(f"nstates {nstates} ;")
        lines.append(f"xread '{title}'")
        lines.append(f"{nchar_total} {ntax}")

        for b in self.blocks:
            lines.append(self._block_tag(b))
            for tax in self.master_taxa:
                row = self._row_for_taxon(block=b, taxon=tax)
                lines.append(f"{tax}\t{row}")

        lines.append(";")
        lines.append("")
        lines.append("ccode-.;")
        lines.append("cname")
        for i, b in enumerate(self.blocks, start=1):
            safe_block_name = sanitize_taxon(b.name)
            lines.append(f"[{i} {safe_block_name} ;")
        lines.append(";")
        lines.append("proc/;")

        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))

    def _block_tag(self, b: DataBlock) -> str:
        if b.dtype == "cont":
            return "&[cont]"
        if b.dtype == "num":
            return "&[num]"
        if b.dtype == "dna":
            return "&[dna]"
        if b.dtype == "landmark2d":
            return "&[landmark 2D]"
        return "&[num]"

    def _row_for_taxon(self, block: DataBlock, taxon: str) -> str:
        """
        Return formatted row data (string after taxon).
        Missing taxa are filled appropriately.
        """
        if block.dtype == "dna":
            if taxon in block.rows_by_taxon:
                seq = block.rows_by_taxon[taxon]
                return re.sub(r"\s+", "", seq)
            return block.missing_symbol * block.nchar

        if block.dtype == "landmark2d":
            k = int(block.meta.get("n_landmarks", 0))
            if k <= 0:
                raise ValueError("Landmark block missing n_landmarks metadata.")
            if taxon not in block.rows_by_taxon:
                # missing taxon => 2*k missing values
                return " ".join([block.missing_symbol] * k)

            coords: List[Tuple[float, float]] = block.rows_by_taxon[taxon]
            if len(coords) != k:
                raise ValueError(f"Landmark count mismatch for {taxon}: expected {k}, got {len(coords)}")

            out_vals: List[str] = []
            pairs: List[str] = []
            for (x, y) in coords:
                pairs.append(f"{x:.{self.decimals}f},{y:.{self.decimals}f}")
            return " ".join(pairs)

        # cont/num: stored as list of tokens
        if taxon in block.rows_by_taxon:
            tokens = block.rows_by_taxon[taxon]
        else:
            tokens = [block.missing_symbol] * block.nchar

        if len(tokens) != block.nchar:
            if len(tokens) < block.nchar:
                tokens = tokens + [block.missing_symbol] * (block.nchar - len(tokens))
            else:
                tokens = tokens[: block.nchar]

        if block.dtype == "cont":
            out_tokens: List[str] = []
            for t in tokens:
                if is_missing(t):
                    out_tokens.append(block.missing_symbol)
                else:
                    fv = parse_float_safe(t)
                    if fv is None or math.isnan(fv):
                        out_tokens.append(block.missing_symbol)
                    else:
                        out_tokens.append(f"{fv:.{self.decimals}f}")
            return " ".join(out_tokens)

        out_tokens = [(block.missing_symbol if is_missing(t) else t) for t in tokens]
        return " ".join(out_tokens)


# -------------------------
# GUI
# -------------------------
class TNTBuilderApp(tk.Tk if not HAS_DND else TkinterDnD.Tk):  
    def __init__(self):
        super().__init__()
        self.title("TNT Matrix Builder Liria & Soto-Vivas (2026)")
        self.geometry("1020x590")
        self.minsize(980, 520)

        self.project = Project()

        self._build_ui()
        self._refresh_ui()

    def _build_ui(self):
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        main = ttk.Frame(self, padding=10)
        main.grid(row=0, column=0, sticky="nsew")
        main.columnconfigure(0, weight=3)
        main.columnconfigure(1, weight=2)
        main.rowconfigure(1, weight=1)

        # --- Drop zone / Add files ---
        drop_frame = ttk.Frame(main)
        drop_frame.grid(row=0, column=0, columnspan=2, sticky="ew")
        drop_frame.columnconfigure(0, weight=1)

        self.drop_label = ttk.Label(
            drop_frame,
            text="Drag & drop files here (FASTA / TXT/TSV/CSV / TPS)."
                 + ("" if HAS_DND else "  (Tip: pip install tkinterdnd2 for drag&drop)"),
            relief="ridge",
            padding=12
        )
        self.drop_label.grid(row=0, column=0, sticky="ew", padx=(0, 8), pady=(0, 10))

        if HAS_DND:
            self.drop_label.drop_target_register(DND_FILES)  
            self.drop_label.dnd_bind("<<Drop>>", self._on_drop)  

        btns = ttk.Frame(drop_frame)
        btns.grid(row=0, column=1, sticky="e", pady=(0, 10))
        ttk.Button(btns, text="Add files…", command=self._add_files).grid(row=0, column=0, padx=4)
        ttk.Button(btns, text="Clear", command=self._clear_all).grid(row=0, column=1, padx=4)

        # --- Left: Blocks table ---
        left = ttk.Frame(main)
        left.grid(row=1, column=0, sticky="nsew", padx=(0, 10))
        left.rowconfigure(1, weight=1)
        left.columnconfigure(0, weight=1)

        ttk.Label(left, text="Blocks (in export order):").grid(row=0, column=0, sticky="w", pady=(0, 6))

        self.tree = ttk.Treeview(left, columns=("type", "ntax", "nchar", "status"), show="headings", height=15)
        self.tree.heading("type", text="Type")
        self.tree.heading("ntax", text="#Taxa")
        self.tree.heading("nchar", text="#Chars")
        self.tree.heading("status", text="Status")
        self.tree.column("type", width=100, anchor="center")
        self.tree.column("ntax", width=85, anchor="center")
        self.tree.column("nchar", width=90, anchor="center")
        self.tree.column("status", width=330, anchor="w")
        self.tree.grid(row=1, column=0, sticky="nsew")

        sb = ttk.Scrollbar(left, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        sb.grid(row=1, column=1, sticky="ns")

        ctrl = ttk.Frame(left)
        ctrl.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        ttk.Button(ctrl, text="Remove selected", command=self._remove_selected).grid(row=0, column=0, padx=4)
        ttk.Button(ctrl, text="Move up", command=lambda: self._move_selected(-1)).grid(row=0, column=1, padx=4)
        ttk.Button(ctrl, text="Move down", command=lambda: self._move_selected(1)).grid(row=0, column=2, padx=4)
        ttk.Button(ctrl, text="Rename block…", command=self._rename_selected).grid(row=0, column=3, padx=4)

        # --- Right: Options + Summary ---
        right = ttk.Frame(main)
        right.grid(row=1, column=1, sticky="nsew")
        right.columnconfigure(0, weight=1)

        taxa_box = ttk.LabelFrame(right, text="Taxa merge", padding=10)
        taxa_box.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        self.merge_mode_var = tk.StringVar(value="union")
        ttk.Radiobutton(taxa_box, text="Union (keep all taxa; fill missing with ?)",
                        variable=self.merge_mode_var, value="union",
                        command=self._on_option_change).grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(taxa_box, text="Intersection (only taxa present in all blocks)",
                        variable=self.merge_mode_var, value="intersection",
                        command=self._on_option_change).grid(row=1, column=0, sticky="w")

        self.normalize_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(taxa_box, text="Normalize taxon names (spaces → _)",
                        variable=self.normalize_var,
                        command=self._on_option_change).grid(row=2, column=0, sticky="w", pady=(6, 0))

        ttk.Button(taxa_box, text="Show taxa report…", command=self._show_report).grid(row=3, column=0, sticky="w", pady=(8, 0))

        export_box = ttk.LabelFrame(right, text="Export options", padding=10)
        export_box.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        export_box.columnconfigure(1, weight=1)

        ttk.Label(export_box, text="nstates:").grid(row=0, column=0, sticky="w")
        self.nstates_var = tk.StringVar(value="auto")
        self.nstates_combo = ttk.Combobox(export_box, textvariable=self.nstates_var,
                                          values=["auto", "8", "16", "32"], width=8, state="readonly")
        self.nstates_combo.grid(row=0, column=1, sticky="w")
        self.nstates_combo.bind("<<ComboboxSelected>>", lambda e: self._on_option_change())

        ttk.Label(export_box, text="Decimals (cont & LM):").grid(row=1, column=0, sticky="w", pady=(6, 0))
        self.decimals_var = tk.StringVar(value="6")
        dec_entry = ttk.Entry(export_box, textvariable=self.decimals_var, width=8)
        dec_entry.grid(row=1, column=1, sticky="w", pady=(6, 0))
        dec_entry.bind("<FocusOut>", lambda e: self._on_option_change())

        summary_box = ttk.LabelFrame(right, text="Summary", padding=10)
        summary_box.grid(row=2, column=0, sticky="ew")
        self.summary_lbl = ttk.Label(summary_box, text="", justify="left")
        self.summary_lbl.grid(row=0, column=0, sticky="w")

        ttk.Button(right, text="Export .tnt…", command=self._export).grid(row=3, column=0, sticky="ew", pady=(10, 0))

    # -------------------------
    # GUI actions
    # -------------------------
    def _on_drop(self, event):
        paths = split_dnd_files(event.data)
        self._load_files(paths)

    def _add_files(self):
        paths = filedialog.askopenfilenames(
            title="Select input files",
            filetypes=[
                ("Supported", "*.fasta *.fa *.fas *.txt *.tsv *.csv *.tps"),
                ("FASTA", "*.fasta *.fa *.fas"),
                ("Tables", "*.txt *.tsv *.csv"),
                ("TPS", "*.tps"),
                ("All files", "*.*"),
            ],
        )
        if paths:
            self._load_files(list(paths))

    def _clear_all(self):
        self.project.blocks.clear()
        self.project.rebuild()
        self._refresh_ui()

    def _remove_selected(self):
        sel = self.tree.selection()
        if not sel:
            return
        idx = int(sel[0])
        self.project.remove_block(idx)
        self._refresh_ui()

    def _move_selected(self, direction: int):
        sel = self.tree.selection()
        if not sel:
            return
        idx = int(sel[0])
        self.project.move_block(idx, direction)
        self._refresh_ui()
        new_idx = max(0, min(len(self.project.blocks) - 1, idx + direction))
        self.tree.selection_set(str(new_idx))

    def _rename_selected(self):
        sel = self.tree.selection()
        if not sel:
            return
        idx = int(sel[0])
        b = self.project.blocks[idx]

        dialog = tk.Toplevel(self)
        dialog.title("Rename block")
        dialog.resizable(False, False)
        dialog.transient(self)
        dialog.grab_set()

        ttk.Label(dialog, text="Block name (used in cname):").grid(row=0, column=0, padx=10, pady=(10, 4), sticky="w")
        var = tk.StringVar(value=b.name)
        ent = ttk.Entry(dialog, textvariable=var, width=45)
        ent.grid(row=1, column=0, padx=10, pady=(0, 10))
        ent.focus_set()

        def ok():
            new = var.get().strip()
            if not new:
                messagebox.showwarning("Invalid", "Block name cannot be empty.")
                return
            b.name = new
            self.project.rebuild()
            self._refresh_ui()
            dialog.destroy()

        ttk.Button(dialog, text="OK", command=ok).grid(row=2, column=0, padx=10, pady=(0, 10), sticky="e")

    def _on_option_change(self):
        self.project.union_mode = self.merge_mode_var.get()
        self.project.normalize_names = bool(self.normalize_var.get())

        try:
            dec = int(self.decimals_var.get().strip())
            dec = max(0, min(12, dec))
        except Exception:
            dec = 6
        self.decimals_var.set(str(dec))
        self.project.decimals = dec

        self.project.nstates_mode = self.nstates_var.get()

        self.project.rebuild()
        self._refresh_ui()

    def _show_report(self):
        if not self.project.blocks:
            messagebox.showinfo("Report", "No blocks loaded.")
            return

        self.project.validate()
        rep = self.project.report

        win = tk.Toplevel(self)
        win.title("Taxa report")
        win.geometry("800x500")

        txt = tk.Text(win, wrap="word")
        txt.pack(fill="both", expand=True)
        txt.insert("end", f"Master taxa mode: {self.project.union_mode}\n")
        txt.insert("end", f"Master taxa count: {len(self.project.master_taxa)}\n")
        txt.insert("end", f"Total characters: {self.project.total_nchar()}\n")
        txt.insert("end", "-" * 70 + "\n\n")

        for b in self.project.blocks:
            r = rep.get(b.name, {})
            txt.insert("end", f"BLOCK: {b.name}\n")
            txt.insert("end", f"  Type: {b.dtype} | nchar={b.nchar} | ntax={len(b.taxa)}\n")
            if b.meta:
                txt.insert("end", f"  Meta: {b.meta}\n")
            miss = r.get("missing_taxa", [])
            if miss:
                txt.insert("end", f"  Missing taxa (will be filled): {len(miss)}\n")
                txt.insert("end", "    " + ", ".join(miss[:50]) + (" ..." if len(miss) > 50 else "") + "\n")
            warns = r.get("warnings", [])
            for w in warns:
                txt.insert("end", f"  Warning: {w}\n")
            txt.insert("end", "\n")

        txt.config(state="disabled")

    def _export(self):
        if not self.project.blocks:
            messagebox.showwarning("Export", "Load at least one file first.")
            return
        if not self.project.master_taxa:
            messagebox.showwarning("Export", "Master taxa is empty (try Union mode).")
            return

        out_path = filedialog.asksaveasfilename(
            title="Save TNT file",
            defaultextension=".tnt",
            filetypes=[("TNT files", "*.tnt"), ("All files", "*.*")]
        )
        if not out_path:
            return

        try:
            self.project.export_tnt(out_path=out_path, title="Data saved from TNT")
        except Exception as e:
            messagebox.showerror("Export failed", str(e))
            return

        messagebox.showinfo("Export", f"Saved:\n{out_path}")

    # -------------------------
    # Loading files
    # -------------------------
    def _load_files(self, paths: List[str]):
        errors: List[str] = []
        for p in paths:
            if not p or not os.path.exists(p):
                errors.append(f"File not found: {p}")
                continue
            try:
                block = Parsers.parse(p)
                self.project.add_block(block)
            except Exception as e:
                errors.append(f"{os.path.basename(p)}: {e}")

        self._refresh_ui()

        if errors:
            messagebox.showwarning("Some files could not be loaded", "\n\n".join(errors))

    # -------------------------
    # UI refresh
    # -------------------------
    def _refresh_ui(self):
        for item in self.tree.get_children():
            self.tree.delete(item)

        for idx, b in enumerate(self.project.blocks):
            status = "OK"
            r = self.project.report.get(b.name, {})
            miss = r.get("missing_taxa", [])
            if miss:
                status = f"Missing taxa: {len(miss)} (filled)"
            if b.warnings:
                status = "Warnings: " + "; ".join(b.warnings[:2]) + (" ..." if len(b.warnings) > 2 else "")
            self.tree.insert("", "end", iid=str(idx), values=(b.short_type(), len(b.taxa), b.nchar, status))

        nstates = self.project.choose_nstates() if self.nstates_var.get() == "auto" else int(self.nstates_var.get())
        summary = (
            f"Blocks: {len(self.project.blocks)}\n"
            f"Master taxa: {len(self.project.master_taxa)} ({self.project.union_mode})\n"
            f"Total characters: {self.project.total_nchar()}\n"
            f"nstates: {nstates}\n"
            f"Order: " + " → ".join([f"{b.dtype}" for b in self.project.blocks])
        )
        self.summary_lbl.config(text=summary)
                
        self.project.union_mode = self.merge_mode_var.get()
        self.project.normalize_names = bool(self.normalize_var.get())
        self.project.nstates_mode = self.nstates_var.get()
        try:
            self.project.decimals = int(self.decimals_var.get().strip())
        except Exception:
            self.project.decimals = 6

def main():
    app = TNTBuilderApp()
    app.mainloop()

if __name__ == "__main__":
    main()
