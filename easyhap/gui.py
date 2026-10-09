"""EasyHap Graphical User Interface (GUI).

Provides a modern, tabbed interface organized into three primary workflows:
  1. Haplotype Analysis (easyhap analyze)
  2. Multi-Locus Haplotype Aggregation (easyhap aggregate)
  3. All-Sample Haplotype Heatmap (easyhap heatmap)
"""

from __future__ import annotations

import os
import subprocess
import sys
import threading
from typing import Optional

try:
    import tkinter as tk
    from tkinter import colorchooser, filedialog, messagebox, ttk
except ImportError:
    tk = None
    ttk = None

from . import __version__
from .core import run_analysis
from .vcf_reader import PlinkBedReader, Plink2PgenReader


def _open_folder_in_explorer(path: str) -> None:
    """Open folder in native file manager across Windows, macOS, and Linux."""
    target = os.path.abspath(path)
    if not os.path.exists(target):
        os.makedirs(target, exist_ok=True)
    try:
        if sys.platform == "win32":
            os.startfile(target)
        elif sys.platform == "darwin":
            subprocess.run(["open", target], check=False)
        else:
            subprocess.run(["xdg-open", target], check=False)
    except Exception:
        pass


class EasyHapGUI(tk.Tk if tk else object):
    """Modern scientific workbench GUI for EasyHap."""

    def __init__(self) -> None:
        if tk is None:
            raise RuntimeError("Tkinter is required to run the graphical interface.")
        super().__init__()
        self._job_running = False
        self.title(
            f"EasyHap v{__version__} — Haplotype Analysis & Visualization Workbench"
        )
        self.geometry("1060x890")
        self.minsize(920, 720)
        self._setup_style()
        self._build_ui()

    def _setup_style(self) -> None:
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass

        # Color palette - Fresh, minimalist, comfortable scientific aesthetic
        self.c_bg = "#F8FAFC"  # Clean soft off-white
        self.c_card = "#FFFFFF"  # Crisp white card
        self.c_header = "#0F2942"  # Elegant deep scientific ocean header
        self.c_primary = "#0284C7"  # Fresh cerulean / sky blue
        self.c_primary_hover = "#0369A1"
        self.c_text = "#1E293B"  # Comfortable slate text
        self.c_subtext = "#64748B"  # Muted subtext
        self.c_border = "#E2E8F0"  # Subtle light border
        self.c_accent = "#0D9488"  # Soft teal accent

        self.configure(bg=self.c_bg)

        style.configure("TFrame", background=self.c_bg)
        style.configure("Card.TFrame", background=self.c_card)
        style.configure("TLabelframe", background=self.c_card, relief="groove")
        style.configure(
            "TLabelframe.Label",
            background=self.c_card,
            foreground=self.c_primary,
            font=("Segoe UI", 10, "bold")
            if sys.platform == "win32"
            else ("Arial", 10, "bold"),
        )

        style.configure("Header.TFrame", background=self.c_header)
        style.configure(
            "HeaderTitle.TLabel",
            background=self.c_header,
            foreground="#FFFFFF",
            font=("Segoe UI", 16, "bold")
            if sys.platform == "win32"
            else ("Arial", 16, "bold"),
        )
        style.configure(
            "HeaderSub.TLabel",
            background=self.c_header,
            foreground="#94A3B8",
            font=("Segoe UI", 9) if sys.platform == "win32" else ("Arial", 9),
        )

        style.configure(
            "TLabel",
            background=self.c_card,
            foreground=self.c_text,
            font=("Segoe UI", 9) if sys.platform == "win32" else ("Arial", 9),
        )
        style.configure(
            "Sub.TLabel",
            background=self.c_card,
            foreground=self.c_subtext,
            font=("Segoe UI", 8) if sys.platform == "win32" else ("Arial", 8),
        )
        style.configure(
            "Req.TLabel",
            background=self.c_card,
            foreground="#DC2626",
            font=("Segoe UI", 9, "bold")
            if sys.platform == "win32"
            else ("Arial", 9, "bold"),
        )

        style.configure(
            "Primary.TButton",
            font=("Segoe UI", 10, "bold")
            if sys.platform == "win32"
            else ("Arial", 10, "bold"),
            padding=(16, 7),
        )
        style.configure(
            "Secondary.TButton",
            font=("Segoe UI", 9) if sys.platform == "win32" else ("Arial", 9),
            padding=(8, 4),
        )
        style.configure("TNotebook", background=self.c_bg, tabmargins=[6, 5, 2, 0])
        style.configure(
            "TNotebook.Tab",
            font=("Segoe UI", 10, "bold")
            if sys.platform == "win32"
            else ("Arial", 10, "bold"),
            padding=(14, 6),
        )

    def _file_browser_row(
        self,
        parent,
        row: int,
        label: str,
        var: tk.StringVar,
        required: bool = False,
        directory: bool = False,
        filetypes: Optional[list] = None,
        hint: Optional[str] = None,
    ) -> None:
        """Create a standard label + entry + browse button row."""
        lbl_text = label + (" *" if required else "")
        lbl = ttk.Label(
            parent, text=lbl_text, style="Req.TLabel" if required else "TLabel"
        )
        lbl.grid(row=row, column=0, sticky="w", padx=(10, 6), pady=4)

        entry = ttk.Entry(parent, textvariable=var)
        entry.grid(row=row, column=1, sticky="ew", padx=6, pady=4)

        def on_browse():
            if directory:
                p = filedialog.askdirectory(title=f"Select {label}")
            else:
                p = filedialog.askopenfilename(
                    title=f"Select {label}",
                    filetypes=filetypes or [("All files", "*.*")],
                )
            if p:
                var.set(p)

        btn = ttk.Button(
            parent, text="Browse…", style="Secondary.TButton", command=on_browse
        )
        btn.grid(row=row, column=2, padx=(6, 10), pady=4)

        if hint:
            hint_lbl = ttk.Label(parent, text=hint, style="Sub.TLabel")
            hint_lbl.grid(row=row + 1, column=1, sticky="w", padx=6, pady=(0, 4))

    def _color_row(self, parent, row: int, label: str, var: tk.StringVar) -> None:
        """Create a color swatch and picker row."""
        ttk.Label(parent, text=label).grid(
            row=row, column=0, sticky="w", padx=(10, 6), pady=3
        )
        ent = ttk.Entry(parent, textvariable=var, width=12)
        ent.grid(row=row, column=1, sticky="w", padx=6, pady=3)

        preview = tk.Label(
            parent, text="   ", bg=var.get() or "#FFFFFF", relief="solid", bd=1
        )
        preview.grid(row=row, column=2, sticky="w", padx=4, pady=3)

        def on_pick():
            c = colorchooser.askcolor(color=var.get(), title=f"Choose {label}")[1]
            if c:
                var.set(c)
                preview.configure(bg=c)

        ttk.Button(
            parent, text="Pick…", style="Secondary.TButton", command=on_pick
        ).grid(row=row, column=3, sticky="w", padx=(4, 10), pady=3)

    def _build_ui(self) -> None:
        # Top banner
        header = ttk.Frame(self, style="Header.TFrame", padding=(18, 12))
        header.pack(fill="x")
        ttk.Label(
            header,
            text="EasyHap  Haplotype Analysis & Visualization",
            style="HeaderTitle.TLabel",
        ).pack(anchor="w")
        ttk.Label(
            header,
            text=f"Version {__version__}  |  Multi-Ploidy Regional Analysis, Multi-Locus Aggregation, All-Sample Heatmap & LD Network",
            style="HeaderSub.TLabel",
        ).pack(anchor="w", pady=(2, 0))

        # Main notebook
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True, padx=12, pady=(10, 6))

        tab_ana = ttk.Frame(self.notebook)
        tab_agg = ttk.Frame(self.notebook)
        tab_hm = ttk.Frame(self.notebook)
        tab_help = ttk.Frame(self.notebook)

        self.notebook.add(tab_ana, text=" 🧬 1. Haplotype Analysis ")
        self.notebook.add(tab_agg, text=" 📊 2. Multi-Locus Aggregation ")
        self.notebook.add(tab_hm, text=" 🗺️ 3. All-Sample Heatmap ")
        self.notebook.add(tab_help, text=" 📖 4. Quick Help & Docs ")

        self._build_analysis_tab(tab_ana)
        self._build_aggregate_tab(tab_agg)
        self._build_heatmap_tab(tab_hm)
        self._build_help_tab(tab_help)

        # Bottom persistent console & status
        self._build_bottom_panel()

    def _build_analysis_tab(self, parent: ttk.Frame) -> None:
        """Tab 1: easyhap analyze interface."""
        canvas = tk.Canvas(parent, bg=self.c_bg, highlightthickness=0)
        vscroll = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=vscroll.set)

        canvas.pack(side="left", fill="both", expand=True)
        vscroll.pack(side="right", fill="y")

        # Section 1: Genotype Data Source
        src_box = ttk.LabelFrame(
            scrollable_frame, text=" Genotype Data Source ", padding=10
        )
        src_box.pack(fill="x", padx=12, pady=(8, 6))
        src_box.columnconfigure(1, weight=1)

        self.ana_source = tk.StringVar(value="vcf")
        self.ana_vcf = tk.StringVar()
        self.ana_bfile = tk.StringVar()
        self.ana_pfile = tk.StringVar()
        self.ana_hmp = tk.StringVar()

        r_bar = ttk.Frame(src_box, style="Card.TFrame")
        r_bar.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8), padx=10)

        for val, txt in [
            ("vcf", "VCF / VCF.GZ / BCF"),
            ("bfile", "PLINK 1 BED/BIM/FAM"),
            ("pfile", "PLINK 2 PGEN/PVAR/PSAM"),
            ("hmp", "HapMap (.hmp / .gz)"),
        ]:
            ttk.Radiobutton(
                r_bar,
                text=txt,
                variable=self.ana_source,
                value=val,
                command=self._on_source_switch,
            ).pack(side="left", padx=(0, 16))

        # VCF row
        ttk.Label(src_box, text="VCF / BCF file").grid(
            row=1, column=0, sticky="w", padx=10, pady=4
        )
        self.vcf_entry = ttk.Entry(src_box, textvariable=self.ana_vcf)
        self.vcf_entry.grid(row=1, column=1, sticky="ew", padx=6, pady=4)
        self.vcf_btn = ttk.Button(
            src_box,
            text="Browse…",
            style="Secondary.TButton",
            command=lambda: self._choose_file(
                self.ana_vcf, [("VCF/BCF", "*.vcf *.vcf.gz *.bcf"), ("All", "*.*")]
            ),
        )
        self.vcf_btn.grid(row=1, column=2, padx=10, pady=4)

        # PLINK 1 row
        ttk.Label(src_box, text="PLINK 1 prefix (.bed/.bim/.fam)").grid(
            row=2, column=0, sticky="w", padx=10, pady=4
        )
        self.bfile_entry = ttk.Entry(src_box, textvariable=self.ana_bfile)
        self.bfile_entry.grid(row=2, column=1, sticky="ew", padx=6, pady=4)
        self.bfile_btn = ttk.Button(
            src_box,
            text="Browse trio…",
            style="Secondary.TButton",
            command=self._choose_plink1,
        )
        self.bfile_btn.grid(row=2, column=2, padx=10, pady=4)

        # PLINK 2 row
        ttk.Label(src_box, text="PLINK 2 prefix (.pgen/.pvar/.psam)").grid(
            row=3, column=0, sticky="w", padx=10, pady=4
        )
        self.pfile_entry = ttk.Entry(src_box, textvariable=self.ana_pfile)
        self.pfile_entry.grid(row=3, column=1, sticky="ew", padx=6, pady=4)
        self.pfile_btn = ttk.Button(
            src_box,
            text="Browse pfile…",
            style="Secondary.TButton",
            command=self._choose_plink2,
        )
        self.pfile_btn.grid(row=3, column=2, padx=10, pady=4)

        # HapMap row
        ttk.Label(src_box, text="HapMap file").grid(
            row=4, column=0, sticky="w", padx=10, pady=4
        )
        self.hmp_entry = ttk.Entry(src_box, textvariable=self.ana_hmp)
        self.hmp_entry.grid(row=4, column=1, sticky="ew", padx=6, pady=4)
        self.hmp_btn = ttk.Button(
            src_box,
            text="Browse…",
            style="Secondary.TButton",
            command=lambda: self._choose_file(
                self.ana_hmp, [("HapMap", "*.hmp *.hmp.txt *.hmp.gz"), ("All", "*.*")]
            ),
        )
        self.hmp_btn.grid(row=4, column=2, padx=10, pady=4)

        # Section 2: Target Region & Annotations
        reg_box = ttk.LabelFrame(
            scrollable_frame, text=" Target Region & Metadata ", padding=10
        )
        reg_box.pack(fill="x", padx=12, pady=6)
        reg_box.columnconfigure(1, weight=1)

        self.ana_region = tk.StringVar()
        self.ana_region_file = tk.StringVar()
        self.ana_variant_list = tk.StringVar()
        self.ana_group = tk.StringVar()
        self.ana_traits = tk.StringVar()
        self.ana_gff = tk.StringVar()
        self.ana_outdir = tk.StringVar(value="EasyHap_results")

        ttk.Label(reg_box, text="Region (single) *", style="Req.TLabel").grid(
            row=0, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(reg_box, textvariable=self.ana_region).grid(
            row=0, column=1, sticky="ew", padx=6, pady=4
        )
        ttk.Label(reg_box, text="e.g. Chr4:4072552-4076934", style="Sub.TLabel").grid(
            row=0, column=2, sticky="w", padx=10
        )

        self._file_browser_row(
            reg_box,
            1,
            "Region file (batch)",
            self.ana_region_file,
            hint="One region per line for batch execution (alternative to single region).",
        )
        self._file_browser_row(
            reg_box,
            3,
            "Variant list (-l / --variant-list)",
            self.ana_variant_list,
            hint="File with one variant ID/POS per line, or comma-separated variant IDs.",
        )
        self._file_browser_row(
            reg_box,
            5,
            "Sample group table",
            self.ana_group,
            hint="Tab-delimited Accession and Type (optional latitude/longitude for map pies).",
        )
        self._file_browser_row(
            reg_box,
            7,
            "Phenotypic trait table",
            self.ana_traits,
            hint="Tab-delimited Accession and trait columns.",
        )
        self._file_browser_row(
            reg_box,
            9,
            "Gene annotation (GFF3/GTF)",
            self.ana_gff,
            hint="GFF3/GTF annotation file for gene structure rendering.",
        )
        self._file_browser_row(
            reg_box, 11, "Output directory", self.ana_outdir, directory=True
        )

        # Section 3: Analysis & Quality Control Parameters
        param_box = ttk.LabelFrame(
            scrollable_frame, text=" Analysis & Quality Control ", padding=10
        )
        param_box.pack(fill="x", padx=12, pady=6)
        param_box.columnconfigure(1, weight=1)
        param_box.columnconfigure(3, weight=1)

        self.ana_mode = tk.StringVar(value="genotype")
        self.ana_hetero = tk.StringVar(value="slash")
        self.ana_max_missing = tk.StringVar(value="")
        self.ana_min_maf = tk.StringVar(value="")
        self.ana_min_var = tk.StringVar(value="2")
        self.ana_cluster = tk.StringVar(value="0.15")
        self.ana_gene_feature = tk.StringVar(value="all")
        self.ana_threads = tk.StringVar(value="1")
        self.ana_overlap_flank = tk.StringVar(value="0")

        self.ana_fisher_groups = tk.StringVar()
        self.ana_fisher_alpha = tk.StringVar()
        self.ana_fisher_adjust = tk.StringVar(value="none")

        self.ana_ld = tk.BooleanVar(value=True)
        self.ana_processed = tk.BooleanVar(value=True)
        self.ana_network = tk.BooleanVar(value=True)

        ttk.Label(param_box, text="Analysis mode").grid(
            row=0, column=0, sticky="w", padx=10, pady=4
        )
        self.mode_combo = ttk.Combobox(
            param_box,
            textvariable=self.ana_mode,
            values=["genotype", "copy"],
            state="readonly",
            width=14,
        )
        self.mode_combo.grid(row=0, column=1, sticky="w", padx=6, pady=4)

        ttk.Label(param_box, text="Heterozygous encoding").grid(
            row=0, column=2, sticky="w", padx=10, pady=4
        )
        ttk.Combobox(
            param_box,
            textvariable=self.ana_hetero,
            values=["slash", "iupac", "missing"],
            state="readonly",
            width=14,
        ).grid(row=0, column=3, sticky="w", padx=6, pady=4)

        ttk.Label(param_box, text="Max missing rate (--missing)").grid(
            row=1, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(param_box, textvariable=self.ana_max_missing, width=14).grid(
            row=1, column=1, sticky="w", padx=6, pady=4
        )

        ttk.Label(param_box, text="Min MAF (--maf)").grid(
            row=1, column=2, sticky="w", padx=10, pady=4
        )
        ttk.Entry(param_box, textvariable=self.ana_min_maf, width=14).grid(
            row=1, column=3, sticky="w", padx=6, pady=4
        )

        ttk.Label(param_box, text="Min variants per region").grid(
            row=2, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(param_box, textvariable=self.ana_min_var, width=14).grid(
            row=2, column=1, sticky="w", padx=6, pady=4
        )

        ttk.Label(param_box, text="Cluster threshold").grid(
            row=2, column=2, sticky="w", padx=10, pady=4
        )
        ttk.Entry(param_box, textvariable=self.ana_cluster, width=14).grid(
            row=2, column=3, sticky="w", padx=6, pady=4
        )

        ttk.Label(param_box, text="Fisher groups (e.g. Cultivar,Landrace)").grid(
            row=3, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(param_box, textvariable=self.ana_fisher_groups, width=18).grid(
            row=3, column=1, sticky="w", padx=6, pady=4
        )

        ttk.Label(param_box, text="Fisher alpha / adjustment").grid(
            row=3, column=2, sticky="w", padx=10, pady=4
        )
        f_sub = ttk.Frame(param_box, style="Card.TFrame")
        f_sub.grid(row=3, column=3, sticky="w", padx=6, pady=4)
        ttk.Entry(f_sub, textvariable=self.ana_fisher_alpha, width=7).pack(
            side="left", padx=(0, 4)
        )
        ttk.Combobox(
            f_sub,
            textvariable=self.ana_fisher_adjust,
            values=["none", "bh"],
            state="readonly",
            width=7,
        ).pack(side="left")

        ttk.Label(param_box, text="Gene feature filter").grid(
            row=4, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Combobox(
            param_box,
            textvariable=self.ana_gene_feature,
            values=["all", "exon", "intron", "cds", "utr"],
            state="readonly",
            width=14,
        ).grid(row=4, column=1, sticky="w", padx=6, pady=4)

        ttk.Label(param_box, text="Threads (-t / --threads)").grid(
            row=4, column=2, sticky="w", padx=10, pady=4
        )
        ttk.Entry(param_box, textvariable=self.ana_threads, width=14).grid(
            row=4, column=3, sticky="w", padx=6, pady=4
        )

        ttk.Label(param_box, text="Overlap flank on each side (bp)").grid(
            row=5, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(param_box, textvariable=self.ana_overlap_flank, width=14).grid(
            row=5, column=1, sticky="w", padx=6, pady=4
        )
        ttk.Label(param_box, text="0 = supplied interval only").grid(
            row=5, column=2, columnspan=2, sticky="w", padx=10, pady=4
        )
        chk_sub = ttk.Frame(param_box, style="Card.TFrame")
        chk_sub.grid(row=6, column=0, columnspan=4, sticky="w", padx=10, pady=(6, 2))
        ttk.Checkbutton(
            chk_sub, text="Calculate LD r² & heatmap", variable=self.ana_ld
        ).pack(side="left", padx=(0, 18))
        ttk.Checkbutton(
            chk_sub, text="Generate Haplotype network", variable=self.ana_network
        ).pack(side="left", padx=(0, 18))
        ttk.Checkbutton(
            chk_sub, text="Save processed tokens tables", variable=self.ana_processed
        ).pack(side="left")

        # Section 4: Visualization & Aesthetics
        vis_box = ttk.LabelFrame(
            scrollable_frame, text=" Publication Visualization & Styling ", padding=10
        )
        vis_box.pack(fill="x", padx=12, pady=6)
        vis_box.columnconfigure(1, weight=1)

        self.ana_plot = tk.BooleanVar(value=True)
        self.ana_plot_format = tk.StringVar(value="pdf,png")
        self.ana_level = tk.StringVar(value="hap")
        self.ana_plot_min = tk.StringVar(value="1")
        self.ana_cell_text = tk.StringVar(value="auto")
        self.ana_trait_cols = tk.StringVar()
        self.ana_hap_palette = tk.StringVar(
            value="#4472C4,#70AD47,#ED7D31,#8064A2,#C00000,#4E79A7,#F28E2B"
        )

        self.ana_ref_color = tk.StringVar(value="#70AD47")
        self.ana_alt_color = tk.StringVar(value="#4472C4")
        self.ana_missing_color = tk.StringVar(value="#D9D9D9")
        self.ana_absence_color = tk.StringVar(value="#FFFFFF")

        ttk.Checkbutton(
            vis_box, text="Enable figure generation", variable=self.ana_plot
        ).grid(row=0, column=0, sticky="w", padx=10, pady=4)
        ttk.Label(vis_box, text="Formats (e.g. pdf,png,svg)").grid(
            row=0, column=1, sticky="w", padx=6, pady=4
        )
        ttk.Entry(vis_box, textvariable=self.ana_plot_format, width=16).grid(
            row=0, column=2, sticky="w", padx=10, pady=4
        )

        ttk.Label(vis_box, text="Haplotype level / Min count").grid(
            row=1, column=0, sticky="w", padx=10, pady=4
        )
        lvl_sub = ttk.Frame(vis_box, style="Card.TFrame")
        lvl_sub.grid(row=1, column=1, columnspan=2, sticky="w", padx=6, pady=4)
        ttk.Combobox(
            lvl_sub,
            textvariable=self.ana_level,
            values=["hap", "cluster"],
            state="readonly",
            width=10,
        ).pack(side="left", padx=(0, 8))
        ttk.Label(lvl_sub, text="Min count (plots/tests/candidates):").pack(side="left", padx=(0, 4))
        ttk.Entry(lvl_sub, textvariable=self.ana_plot_min, width=8).pack(side="left")

        self._color_row(vis_box, 2, "Heatmap REF color", self.ana_ref_color)
        self._color_row(vis_box, 3, "Heatmap ALT color", self.ana_alt_color)
        self._color_row(vis_box, 4, "Missing GT color", self.ana_missing_color)
        self._color_row(vis_box, 5, "Absence color", self.ana_absence_color)

        ttk.Label(vis_box, text="Haplotype palette").grid(
            row=6, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(vis_box, textvariable=self.ana_hap_palette).grid(
            row=6, column=1, sticky="ew", padx=6, pady=4
        )
        ttk.Label(vis_box, text="Comma-separated colors", style="Sub.TLabel").grid(
            row=6, column=2, sticky="w", padx=10
        )

        # Run button for Tab 1
        btn_frame = ttk.Frame(scrollable_frame, padding=12)
        btn_frame.pack(fill="x", padx=12, pady=10)
        ttk.Button(
            btn_frame,
            text="🚀  Run Regional Analysis",
            style="Primary.TButton",
            command=self._run_analysis_action,
        ).pack(side="left")

        self._on_source_switch()

    def _build_aggregate_tab(self, parent: ttk.Frame) -> None:
        """Tab 2: easyhap aggregate interface."""
        canvas = tk.Canvas(parent, bg=self.c_bg, highlightthickness=0)
        vscroll = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=vscroll.set)

        canvas.pack(side="left", fill="both", expand=True)
        vscroll.pack(side="right", fill="y")

        # Intro card
        intro = ttk.LabelFrame(
            scrollable_frame,
            text=" Multi-Locus Haplotype Aggregation Overview ",
            padding=12,
        )
        intro.pack(fill="x", padx=12, pady=(8, 6))
        ttk.Label(
            intro,
            text="Aggregate haplotypes from 2 or more regions produced by 'easyhap analyze', identify all multi-locus\n"
            "haplotype combinations per accession, perform ANOVA, Kruskal-Wallis & pairwise statistical tests with\n"
            "Benjamini-Hochberg FDR correction, rank optimal combinations first, and plot publication boxplots.",
            style="Sub.TLabel",
        ).pack(anchor="w")

        # Section 1: Source & Regions
        sec1 = ttk.LabelFrame(
            scrollable_frame,
            text=" Analysis Directory & Regions to Aggregate ",
            padding=10,
        )
        sec1.pack(fill="x", padx=12, pady=6)
        sec1.columnconfigure(1, weight=1)

        self.agg_indir = tk.StringVar(value="EasyHap_results")
        self.agg_regions = tk.StringVar()
        self.agg_region_file = tk.StringVar()
        self.agg_outdir = tk.StringVar(value="EasyHap_aggregate")
        self.agg_group = tk.StringVar()
        self.agg_map_style = tk.StringVar(value="auto")

        self._file_browser_row(
            sec1,
            0,
            "Analyze results directory *",
            self.agg_indir,
            required=True,
            directory=True,
            hint="Directory containing completed easyhap analyze results (with .HapGroup.tsv files).",
        )

        ttk.Label(sec1, text="Regions (comma-separated)").grid(
            row=2, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(sec1, textvariable=self.agg_regions).grid(
            row=2, column=1, sticky="ew", padx=6, pady=4
        )
        ttk.Label(
            sec1, text="e.g. Chr1:100-500, Chr2:200-600", style="Sub.TLabel"
        ).grid(row=2, column=2, sticky="w", padx=10)

        self._file_browser_row(
            sec1,
            3,
            "Region combination file (--region-file)",
            self.agg_region_file,
            hint="Batch file: each line is one combination, with 2 or more regions separated by TAB.",
        )
        self._file_browser_row(
            sec1, 5, "Aggregate output directory", self.agg_outdir, directory=True
        )
        self._file_browser_row(
            sec1, 7, "Sample groups / coordinates (optional)", self.agg_group,
            hint="Sample, group, latitude, longitude, optional location; coordinates enable combination maps.",
        )
        ttk.Label(sec1, text="Combination map style").grid(row=9, column=0, sticky="w", padx=10, pady=4)
        ttk.Combobox(
            sec1, textvariable=self.agg_map_style, values=["auto", "pie", "bar", "both", "none"],
            state="readonly", width=12,
        ).grid(row=9, column=1, sticky="w", padx=6, pady=4)

        # Section 2: Trait Association & Plotting
        sec2 = ttk.LabelFrame(
            scrollable_frame,
            text=" Phenotypic Trait Association & Boxplots ",
            padding=10,
        )
        sec2.pack(fill="x", padx=12, pady=6)
        sec2.columnconfigure(1, weight=1)

        self.agg_traits = tk.StringVar()
        self.agg_trait_cols = tk.StringVar()
        self.agg_min_count = tk.StringVar(value="1")
        self.agg_plot = tk.BooleanVar(value=True)
        self.agg_plot_format = tk.StringVar(value="pdf,png")
        self.agg_palette = tk.StringVar(
            value="#4472C4,#70AD47,#ED7D31,#8064A2,#C00000,#4E79A7,#F28E2B"
        )

        self._file_browser_row(
            sec2,
            0,
            "External trait table (optional)",
            self.agg_traits,
            hint="Optional TAB-delimited trait file (if omitted, traits are read from analysis HapGroup.tsv results).",
        )

        ttk.Label(sec2, text="Trait columns to analyze").grid(
            row=2, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(sec2, textvariable=self.agg_trait_cols).grid(
            row=2, column=1, sticky="ew", padx=6, pady=4
        )
        ttk.Label(
            sec2,
            text="comma-separated; blank = all available traits",
            style="Sub.TLabel",
        ).grid(row=2, column=2, sticky="w", padx=10)

        ttk.Label(sec2, text="Minimum count for plots").grid(
            row=3, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(sec2, textvariable=self.agg_min_count, width=12).grid(
            row=3, column=1, sticky="w", padx=6, pady=4
        )

        ttk.Checkbutton(
            sec2, text="Generate combination figures and maps", variable=self.agg_plot
        ).grid(row=4, column=0, sticky="w", padx=10, pady=4)
        ttk.Label(sec2, text="Figure formats (pdf,png,svg)").grid(
            row=4, column=1, sticky="w", padx=6, pady=4
        )
        ttk.Entry(sec2, textvariable=self.agg_plot_format, width=16).grid(
            row=4, column=2, sticky="w", padx=10, pady=4
        )

        ttk.Label(sec2, text="Color palette for combinations").grid(
            row=5, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(sec2, textvariable=self.agg_palette).grid(
            row=5, column=1, sticky="ew", padx=6, pady=4
        )
        ttk.Label(sec2, text="comma-separated hex codes", style="Sub.TLabel").grid(
            row=5, column=2, sticky="w", padx=10
        )

        self.agg_test_min = tk.StringVar(value="2")
        self.agg_candidate_min = tk.StringVar(value="5")
        self.agg_direction = tk.StringVar(value="higher")
        self.agg_fdr_scope = tk.StringVar(value="run")
        self.agg_phase_link = tk.BooleanVar(value=False)
        self.agg_interaction = tk.BooleanVar(value=False)
        for row, label, var in [
            (6, "Minimum phenotype N for tests", self.agg_test_min),
            (7, "Minimum phenotype N for candidates", self.agg_candidate_min),
        ]:
            ttk.Label(sec2, text=label).grid(
                row=row, column=0, sticky="w", padx=10, pady=4
            )
            ttk.Entry(sec2, textvariable=var, width=12).grid(
                row=row, column=1, sticky="w", padx=6, pady=4
            )
        ttk.Label(sec2, text="Trait ranking direction").grid(
            row=8, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Combobox(
            sec2,
            textvariable=self.agg_direction,
            values=["higher", "lower"],
            state="readonly",
            width=12,
        ).grid(row=8, column=1, sticky="w", padx=6)
        ttk.Label(sec2, text="FDR scope").grid(
            row=9, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Combobox(
            sec2,
            textvariable=self.agg_fdr_scope,
            values=["run", "trait"],
            state="readonly",
            width=12,
        ).grid(row=9, column=1, sticky="w", padx=6)
        ttk.Checkbutton(
            sec2,
            text="Connect known phase blocks (same chromosome and PS)",
            variable=self.agg_phase_link,
        ).grid(row=10, column=0, columnspan=3, sticky="w", padx=10, pady=4)
        ttk.Checkbutton(
            sec2,
            text="Two-region categorical interaction model",
            variable=self.agg_interaction,
        ).grid(row=11, column=0, columnspan=3, sticky="w", padx=10, pady=4)
        ttk.Label(
            sec2,
            text="Exploratory comparisons; candidate ranks describe trait means.",
            style="Sub.TLabel",
        ).grid(row=12, column=0, columnspan=3, sticky="w", padx=10, pady=4)

        btn_frame = ttk.Frame(scrollable_frame, padding=12)
        btn_frame.pack(fill="x", padx=12, pady=10)
        ttk.Button(
            btn_frame,
            text="📊  Run Multi-Locus Aggregation",
            style="Primary.TButton",
            command=self._run_aggregate_action,
        ).pack(side="left")

    def _build_heatmap_tab(self, parent: ttk.Frame) -> None:
        """Tab 3: easyhap heatmap interface."""
        canvas = tk.Canvas(parent, bg=self.c_bg, highlightthickness=0)
        vscroll = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=vscroll.set)

        canvas.pack(side="left", fill="both", expand=True)
        vscroll.pack(side="right", fill="y")

        intro = ttk.LabelFrame(
            scrollable_frame,
            text=" All-Sample Haplotype Color Heatmap Overview ",
            padding=12,
        )
        intro.pack(fill="x", padx=12, pady=(8, 6))
        ttk.Label(
            intro,
            text="Generate a clean, high-density all-sample haplotype color heatmap with population group annotation\n"
            "bars and within-group hierarchical sample clustering. Displays pure color cells (REF, ALT, HET, Missing, Absence)\n"
            "without text clutter, optimized for large germplasm collections.",
            style="Sub.TLabel",
        ).pack(anchor="w")

        sec1 = ttk.LabelFrame(
            scrollable_frame, text=" Target Region & Directory ", padding=10
        )
        sec1.pack(fill="x", padx=12, pady=6)
        sec1.columnconfigure(1, weight=1)

        self.hm_indir = tk.StringVar(value="EasyHap_results")
        self.hm_region = tk.StringVar()
        self.hm_region_file = tk.StringVar()
        self.hm_group = tk.StringVar()
        self.hm_outdir = tk.StringVar()

        self._file_browser_row(
            sec1,
            0,
            "Analyze results directory *",
            self.hm_indir,
            required=True,
            directory=True,
        )

        ttk.Label(sec1, text="Region to plot (single)").grid(
            row=2, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(sec1, textvariable=self.hm_region).grid(
            row=2, column=1, sticky="ew", padx=6, pady=4
        )
        ttk.Label(
            sec1,
            text="e.g. Chr4:4072552-4076934 (auto if 1 region)",
            style="Sub.TLabel",
        ).grid(row=2, column=2, sticky="w", padx=10)

        self._file_browser_row(
            sec1,
            3,
            "Region file (--region-file)",
            self.hm_region_file,
            hint="Batch file: one region per line for generating heatmaps across multiple regions.",
        )
        self._file_browser_row(
            sec1,
            5,
            "Sample group table (optional)",
            self.hm_group,
            hint="Optional group file; if omitted, automatically read from HapGroup.tsv.",
        )
        self._file_browser_row(
            sec1,
            7,
            "Output directory",
            self.hm_outdir,
            directory=True,
            hint="Leave blank to save directly into the results directory.",
        )

        sec2 = ttk.LabelFrame(
            scrollable_frame, text=" Clustering & Color Customization ", padding=10
        )
        sec2.pack(fill="x", padx=12, pady=6)
        sec2.columnconfigure(1, weight=1)

        self.hm_cluster_method = tk.StringVar(value="average")
        self.hm_plot_format = tk.StringVar(value="pdf,png")

        self.hm_ref_color = tk.StringVar(value="#70AD47")
        self.hm_alt_color = tk.StringVar(value="#4472C4")
        self.hm_het_color = tk.StringVar(value="#F4B183")
        self.hm_missing_color = tk.StringVar(value="#D9D9D9")
        self.hm_absence_color = tk.StringVar(value="#FFFFFF")
        self.hm_palette = tk.StringVar(
            value="#4472C4,#70AD47,#ED7D31,#8064A2,#C00000,#4E79A7,#F28E2B"
        )

        ttk.Label(sec2, text="Within-group clustering method").grid(
            row=0, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Combobox(
            sec2,
            textvariable=self.hm_cluster_method,
            values=["average", "complete", "single", "ward"],
            state="readonly",
            width=14,
        ).grid(row=0, column=1, sticky="w", padx=6, pady=4)

        ttk.Label(sec2, text="Figure formats (pdf,png,svg)").grid(
            row=1, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(sec2, textvariable=self.hm_plot_format, width=16).grid(
            row=1, column=1, sticky="w", padx=6, pady=4
        )

        self._color_row(sec2, 2, "REF allele color", self.hm_ref_color)
        self._color_row(sec2, 3, "ALT allele color", self.hm_alt_color)
        self._color_row(sec2, 4, "HET (heterozygous) color", self.hm_het_color)
        self._color_row(sec2, 5, "Missing genotype color", self.hm_missing_color)
        self._color_row(
            sec2, 6, "Absence / structural absence color", self.hm_absence_color
        )

        ttk.Label(sec2, text="Group annotation palette").grid(
            row=7, column=0, sticky="w", padx=10, pady=4
        )
        ttk.Entry(sec2, textvariable=self.hm_palette).grid(
            row=7, column=1, sticky="ew", padx=6, pady=4
        )

        btn_frame = ttk.Frame(scrollable_frame, padding=12)
        btn_frame.pack(fill="x", padx=12, pady=10)
        ttk.Button(
            btn_frame,
            text="🎨  Generate All-Sample Heatmap",
            style="Primary.TButton",
            command=self._run_heatmap_action,
        ).pack(side="left")

    def _build_help_tab(self, parent: ttk.Frame) -> None:
        """Tab 4: Documentation and usage hints."""
        txt = tk.Text(
            parent,
            wrap="word",
            font=("Consolas", 10),
            padx=16,
            pady=16,
            bg="#FFFFFF",
            fg=self.c_text,
        )
        txt.pack(fill="both", expand=True)

        help_content = f"""EasyHap {__version__} — Research Edition Quick Reference

1. WORKFLOW MODES
--------------------------------------------------------------------------------
EasyHap supports three dedicated operations:
  • Haplotype Analysis (easyhap analyze):
    Performs regional haplotype construction, clustering, Fisher's tests, LD block
    calculation, MST network modeling, and publication-ready figures.
  • Multi-Locus Haplotype Aggregation (easyhap aggregate):
    Aggregates haplotypes from 2 or more regions, associates multi-locus combinations
    with traits (ANOVA, Kruskal-Wallis, pairwise comparisons with BH correction),
    ranks the most favorable combinations first, and generates boxplots.
  • All-Sample Haplotype Heatmap (easyhap heatmap):
    Produces a pure-color genotype heatmap across all individuals with population
    group annotation bars and within-group hierarchical sample clustering.

2. SUPPORTED GENOTYPE FORMATS (Direct .gz support enabled)
--------------------------------------------------------------------------------
  • VCF / BCF: *.vcf, *.vcf.gz, *.bcf
  • PLINK 1: *.bed, *.bim, *.fam (or common prefix; automatically derived)
  • PLINK 2: *.pgen, *.pvar, *.psam (or common prefix; automatically derived)
  • HapMap: *.hmp, *.hmp.txt, *.hmp.gz

3. INPUT REQUIREMENTS
--------------------------------------------------------------------------------
  • Region string: FORMAT 'Chr:start-end', e.g. 'Chr4:4072552-4076934'.
  • Region file for analyze/heatmap: text file with one region per line.
  • Region file for aggregate: text file where each line is a combination of
    2 or more regions separated by TAB (\\t).
  • Sample group table: Tab-delimited (Accession <TAB> Group [<TAB> Lat <TAB> Lon]).
  • Phenotypic traits: Tab-delimited (Accession <TAB> Trait1 <TAB> Trait2 ...).
  • Gene annotation: GFF3 / GTF matching the reference genome.

4. INDEL & SV TOKEN ENCODING
--------------------------------------------------------------------------------
  • Unique-length indels at a site are encoded directly as 'DEL xxbp' or 'INS xxbp'.
  • Equal-length but sequence-distinct indels are prefixed with 'ALT1:DEL xxbp', etc.
  • Standard single-base and heterozygous SNP boxes remain uniform in width (1.0);
    only columns containing structural variant tokens expand to fit the text.
"""
        txt.insert("1.0", help_content)
        txt.configure(state="disabled")

    def _build_bottom_panel(self) -> None:
        """Bottom console log and status bar."""
        frame = ttk.Frame(self, padding=(14, 6))
        frame.pack(fill="both")

        status_bar = ttk.Frame(frame)
        status_bar.pack(fill="x", pady=(0, 6))

        self.status = tk.StringVar(value="Ready")
        ttk.Label(status_bar, text="Status:", font=("Arial", 9, "bold")).pack(
            side="left"
        )
        ttk.Label(
            status_bar,
            textvariable=self.status,
            foreground=self.c_primary,
            font=("Arial", 9, "bold"),
        ).pack(side="left", padx=(6, 16))

        ttk.Button(
            status_bar,
            text="📁 Open Output Folder",
            style="Secondary.TButton",
            command=self._open_current_outdir,
        ).pack(side="right")
        ttk.Button(
            status_bar,
            text="🧹 Clear Log",
            style="Secondary.TButton",
            command=self._clear_log,
        ).pack(side="right", padx=(0, 8))

        log_container = ttk.Frame(frame)
        log_container.pack(fill="both", expand=True)

        self.log_txt = tk.Text(
            log_container,
            height=8,
            font=("Consolas", 9) if sys.platform == "win32" else ("Menlo", 9),
            bg="#FFFFFF",
            fg="#1E293B",
            insertbackground="#0284C7",
            relief="solid",
            bd=1,
            highlightthickness=1,
            highlightbackground="#CBD5E1",
            highlightcolor="#0284C7",
            padx=10,
            pady=8,
        )
        self.log_txt.tag_configure("info", foreground="#0284C7")
        self.log_txt.tag_configure(
            "warn",
            foreground="#D97706",
            font=("Consolas", 9, "bold")
            if sys.platform == "win32"
            else ("Menlo", 9, "bold"),
        )
        self.log_txt.tag_configure(
            "error",
            foreground="#DC2626",
            font=("Consolas", 9, "bold")
            if sys.platform == "win32"
            else ("Menlo", 9, "bold"),
        )
        self.log_txt.tag_configure(
            "success",
            foreground="#16A34A",
            font=("Consolas", 9, "bold")
            if sys.platform == "win32"
            else ("Menlo", 9, "bold"),
        )
        self.log_txt.tag_configure("time", foreground="#94A3B8")

        scroll = ttk.Scrollbar(
            log_container, orient="vertical", command=self.log_txt.yview
        )
        self.log_txt.configure(yscrollcommand=scroll.set)

        self.log_txt.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")

    def _on_source_switch(self) -> None:
        src = self.ana_source.get()
        if src in {"bfile", "pfile"}:
            self.mode_combo.configure(values=["genotype"])
            self.ana_mode.set("genotype")
        else:
            self.mode_combo.configure(values=["genotype", "copy"])

        self.vcf_entry.configure(state="normal" if src == "vcf" else "disabled")
        self.vcf_btn.configure(state="normal" if src == "vcf" else "disabled")
        self.bfile_entry.configure(state="normal" if src == "bfile" else "disabled")
        self.bfile_btn.configure(state="normal" if src == "bfile" else "disabled")
        self.pfile_entry.configure(state="normal" if src == "pfile" else "disabled")
        self.pfile_btn.configure(state="normal" if src == "pfile" else "disabled")
        self.hmp_entry.configure(state="normal" if src == "hmp" else "disabled")
        self.hmp_btn.configure(state="normal" if src == "hmp" else "disabled")

    def _choose_file(self, var: tk.StringVar, filetypes: list) -> None:
        p = filedialog.askopenfilename(filetypes=filetypes)
        if p:
            var.set(p)

    def _choose_plink1(self) -> None:
        p = filedialog.askopenfilename(
            title="Select PLINK 1 BED/BIM/FAM file",
            filetypes=[("PLINK 1", "*.bed *.bim *.fam"), ("All", "*.*")],
        )
        if p:
            prefix = PlinkBedReader._normalize_prefix(p)
            self.ana_bfile.set(prefix)
            self.ana_source.set("bfile")
            self._on_source_switch()

    def _choose_plink2(self) -> None:
        p = filedialog.askopenfilename(
            title="Select PLINK 2 PGEN/PVAR/PSAM file",
            filetypes=[
                ("PLINK 2", "*.pgen *.pvar *.psam *.pgen.gz *.pvar.gz *.psam.gz"),
                ("All", "*.*"),
            ],
        )
        if p:
            prefix = Plink2PgenReader._normalize_prefix(p)
            self.ana_pfile.set(prefix)
            self.ana_source.set("pfile")
            self._on_source_switch()

    def _log(self, text: str) -> None:
        from datetime import datetime

        time_str = datetime.now().strftime("[%H:%M:%S] ")
        self.log_txt.insert("end", time_str, "time")
        upper = text.upper()
        if "[ERROR]" in upper or "ERROR" in upper or "FAILED" in upper:
            self.log_txt.insert("end", text + "\n", "error")
        elif "[WARNING]" in upper or "[SKIP]" in upper or "WARN" in upper:
            self.log_txt.insert("end", text + "\n", "warn")
        elif (
            "[INFO]" in upper
            or "[EASYHAP]" in upper
            or "[ANALYSIS]" in upper
            or "[AGGREGATE]" in upper
        ):
            self.log_txt.insert("end", text + "\n", "info")
        elif "SUCCESS" in upper or "COMPLETE" in upper or "FINISHED" in upper:
            self.log_txt.insert("end", text + "\n", "success")
        else:
            self.log_txt.insert("end", text + "\n")
        self.log_txt.see("end")

    def _clear_log(self) -> None:
        self.log_txt.delete("1.0", "end")

    def _open_current_outdir(self) -> None:
        current_tab = self.notebook.index(self.notebook.select())
        if current_tab == 0:
            target = self.ana_outdir.get().strip() or "EasyHap_results"
        elif current_tab == 1:
            target = self.agg_outdir.get().strip() or "EasyHap_aggregate"
        else:
            target = (
                self.hm_outdir.get().strip()
                or self.hm_indir.get().strip()
                or "EasyHap_results"
            )
        _open_folder_in_explorer(target)

    def _start_worker(self, worker):
        self._job_running = True

        def guarded():
            try:
                worker()
            finally:
                self.after(0, lambda: setattr(self, "_job_running", False))

        threading.Thread(target=guarded, daemon=True).start()

    def _run_analysis_action(self) -> None:
        if self._job_running:
            messagebox.showinfo(
                "Analysis Running", "Wait for the current analysis to finish."
            )
            return
        src = self.ana_source.get()
        vcf = bfile = pfile = hmp = None
        if src == "vcf":
            vcf = self.ana_vcf.get().strip()
            if not vcf:
                messagebox.showerror("Missing Input", "Please provide a VCF/BCF file.")
                return
        elif src == "bfile":
            bfile = self.ana_bfile.get().strip()
            if not bfile:
                messagebox.showerror(
                    "Missing Input", "Please provide a PLINK 1 BED/BIM/FAM prefix."
                )
                return
        elif src == "pfile":
            pfile = self.ana_pfile.get().strip()
            if not pfile:
                messagebox.showerror(
                    "Missing Input", "Please provide a PLINK 2 PGEN/PVAR/PSAM prefix."
                )
                return
        elif src == "hmp":
            hmp = self.ana_hmp.get().strip()
            if not hmp:
                messagebox.showerror("Missing Input", "Please provide a HapMap file.")
                return

        reg = self.ana_region.get().strip()
        reg_file = self.ana_region_file.get().strip()
        vl = self.ana_variant_list.get().strip()
        if not vl and bool(reg) == bool(reg_file):
            messagebox.showerror(
                "Region Error",
                "Provide at least one of Region, Region file, or Variant list.",
            )
            return

        outdir = self.ana_outdir.get().strip() or "EasyHap_results"
        try:
            max_missing = (
                float(self.ana_max_missing.get().strip())
                if self.ana_max_missing.get().strip()
                else None
            )
            min_maf = (
                float(self.ana_min_maf.get().strip())
                if self.ana_min_maf.get().strip()
                else None
            )
            overlap_flank_bp = int(self.ana_overlap_flank.get().strip() or "0")
            if overlap_flank_bp < 0:
                raise ValueError("Overlap flank must be a nonnegative integer (bp)")
            if overlap_flank_bp and not (reg or reg_file):
                raise ValueError("Positive overlap flank requires a region or region file")

            formats = [
                f.strip() for f in self.ana_plot_format.get().split(",") if f.strip()
            ] or ["pdf"]

            kwargs = dict(
                vcf_path=vcf,
                bfile=bfile,
                pfile=pfile,
                hmp=hmp,
                group_file=self.ana_group.get().strip() or None,
                outdir=outdir,
                region=reg or None,
                region_file=reg_file or None,
                mode=self.ana_mode.get(),
                hetero_policy=self.ana_hetero.get(),
                trait_file=self.ana_traits.get().strip() or None,
                gff_file=self.ana_gff.get().strip() or None,
                gene_feature=self.ana_gene_feature.get(),
                variant_list=self.ana_variant_list.get().strip() or None,
                max_missing=max_missing,
                min_maf=min_maf,
                overlap_flank_bp=overlap_flank_bp,
                threads=int(self.ana_threads.get().strip() or "1"),
                min_variants=int(self.ana_min_var.get() or "2"),
                cluster_threshold=float(self.ana_cluster.get() or "0.15"),
                make_ld=self.ana_ld.get(),
                make_network=self.ana_network.get(),
                write_processed=self.ana_processed.get(),
                make_plots=self.ana_plot.get(),
                plot_formats=formats,
                plot_hap_level=self.ana_level.get(),
                plot_min_count=int(self.ana_plot_min.get() or "1"),
                ref_color=self.ana_ref_color.get(),
                alt_color=self.ana_alt_color.get(),
                missing_color=self.ana_missing_color.get(),
                absence_color=self.ana_absence_color.get(),
                hap_palette=[
                    c.strip()
                    for c in self.ana_hap_palette.get().split(",")
                    if c.strip()
                ]
                or None,
            )

        except (ValueError, TypeError) as exc:
            messagebox.showerror("Invalid Analysis Parameter", str(exc))
            return

        self.status.set("Running Regional Analysis…")
        self._log(
            f"[Analysis] Starting EasyHap analysis for {reg or reg_file} -> {outdir}"
        )

        def worker():
            try:
                res = run_analysis(**kwargs)
                msg = f"[Analysis] Finished {len(res)} region(s) successfully."
                self.after(0, lambda: self.status.set(f"Analysis {res.status}"))
                self.after(0, lambda: self._log(msg))
                if res.runtime_s is not None:
                    self.after(
                        0,
                        lambda: self._log(
                            f"[Performance] Runtime: {res.runtime_s:.2f}s | Peak RAM: {res.peak_ram_gb:.3f} GB"
                        ),
                    )
                self.after(
                    0,
                    lambda: messagebox.showinfo(
                        "Success",
                        f"Regional analysis status: {res.status}\nResults in: {outdir}",
                    ),
                )
            except Exception as e:
                self.after(0, lambda: self.status.set("Analysis Error"))
                self.after(0, lambda e=e: self._log(f"[Error] {e}"))
                self.after(
                    0, lambda e=e: messagebox.showerror("Execution Error", str(e))
                )

        self._start_worker(worker)

    def _run_aggregate_action(self) -> None:
        if self._job_running:
            messagebox.showinfo(
                "Analysis Running", "Wait for the current analysis to finish."
            )
            return
        indir = self.agg_indir.get().strip()
        if not indir or not os.path.isdir(indir):
            messagebox.showerror(
                "Missing Directory",
                "Please specify a valid analysis results directory.",
            )
            return

        regions_raw = self.agg_regions.get().strip()
        reg_file = self.agg_region_file.get().strip()
        if not regions_raw and not reg_file:
            messagebox.showerror(
                "Missing Regions",
                "Please provide regions (comma-separated) or a region combination file.",
            )
            return

        regions_list = None
        if regions_raw:
            regions_list = [
                r.strip() for r in regions_raw.replace(";", ",").split(",") if r.strip()
            ]
            if len(regions_list) < 2:
                regions_list = [r.strip() for r in regions_raw.split() if r.strip()]

        outdir = self.agg_outdir.get().strip() or "EasyHap_aggregate"
        formats = [
            f.strip() for f in self.agg_plot_format.get().split(",") if f.strip()
        ] or ["pdf"]

        from .aggregate import run_aggregate

        try:
            aggregate_kwargs = dict(
                indir=indir,
                regions=regions_list,
                region_file=reg_file or None,
                outdir=outdir,
                trait_file=self.agg_traits.get().strip() or None,
                group_file=self.agg_group.get().strip() or None,
                map_style=self.agg_map_style.get(),
                trait_cols=[
                    c.strip() for c in self.agg_trait_cols.get().split(",") if c.strip()
                ]
                or None,
                make_plots=self.agg_plot.get(),
                plot_formats=formats,
                min_count=int(self.agg_min_count.get() or "1"),
                test_min_count=int(self.agg_test_min.get() or "2"),
                candidate_min_count=int(self.agg_candidate_min.get() or "5"),
                trait_direction=self.agg_direction.get(),
                fdr_scope=self.agg_fdr_scope.get(),
                phase_link=self.agg_phase_link.get(),
                interaction=self.agg_interaction.get(),
                hap_palette=[
                    c.strip() for c in self.agg_palette.get().split(",") if c.strip()
                ]
                or None,
            )
        except (ValueError, TypeError) as exc:
            messagebox.showerror("Invalid Aggregate Parameter", str(exc))
            return

        self.status.set("Running Aggregation…")
        self._log(
            f"[Aggregate] Aggregating multi-locus combinations from {indir} -> {outdir}"
        )

        def worker():
            try:
                res = run_aggregate(**aggregate_kwargs)
                self.after(0, lambda: self.status.set("Aggregation Complete"))
                self.after(
                    0,
                    lambda: self._log(
                        f"[Aggregate] Aggregation finished successfully. Results in: {outdir}"
                    ),
                )
                self.after(
                    0,
                    lambda: messagebox.showinfo(
                        "Success",
                        f"Multi-locus aggregation completed!\nOutputs saved to: {outdir}",
                    ),
                )
            except Exception as e:
                self.after(0, lambda: self.status.set("Aggregation Error"))
                self.after(0, lambda e=e: self._log(f"[Error] {e}"))
                self.after(
                    0, lambda e=e: messagebox.showerror("Aggregation Error", str(e))
                )

        self._start_worker(worker)

    def _run_heatmap_action(self) -> None:
        if self._job_running:
            messagebox.showinfo(
                "Analysis Running", "Wait for the current analysis to finish."
            )
            return
        indir = self.hm_indir.get().strip()
        if not indir or not os.path.isdir(indir):
            messagebox.showerror(
                "Missing Directory",
                "Please specify a valid analysis results directory.",
            )
            return

        reg = self.hm_region.get().strip() or None
        reg_file = self.hm_region_file.get().strip() or None
        outdir = self.hm_outdir.get().strip() or indir
        formats = [
            f.strip() for f in self.hm_plot_format.get().split(",") if f.strip()
        ] or ["pdf"]

        from .all_sample_heatmap import run_all_sample_heatmap

        heatmap_kwargs = dict(
            indir=indir,
            region=reg,
            region_file=reg_file,
            group_file=self.hm_group.get().strip() or None,
            outdir=outdir,
            plot_formats=formats,
            cluster_method=self.hm_cluster_method.get(),
            ref_color=self.hm_ref_color.get(),
            alt_color=self.hm_alt_color.get(),
            het_color=self.hm_het_color.get(),
            missing_color=self.hm_missing_color.get(),
            absence_color=self.hm_absence_color.get(),
            hap_palette=[
                c.strip() for c in self.hm_palette.get().split(",") if c.strip()
            ]
            or None,
        )

        self.status.set("Generating All-Sample Heatmap…")
        self._log(
            f"[Heatmap] Generating all-sample heatmap for {reg or reg_file or 'auto-detected region'}..."
        )

        def worker():
            try:
                files = run_all_sample_heatmap(**heatmap_kwargs)
                self.after(0, lambda: self.status.set("Heatmap Complete"))
                self.after(
                    0,
                    lambda: self._log(
                        f"[Heatmap] Generated {len(files)} heatmap(s) in: {outdir}"
                    ),
                )
                self.after(
                    0,
                    lambda: messagebox.showinfo(
                        "Success",
                        f"All-sample heatmap generated successfully!\nFiles saved in: {outdir}",
                    ),
                )
            except Exception as e:
                self.after(0, lambda: self.status.set("Heatmap Error"))
                self.after(0, lambda e=e: self._log(f"[Error] {e}"))
                self.after(0, lambda e=e: messagebox.showerror("Heatmap Error", str(e)))

        self._start_worker(worker)


def main() -> None:
    if tk is None:
        print(
            "[Error] Tkinter is not installed or available in this Python environment.",
            file=sys.stderr,
        )
        print(
            "To run the GUI, please run on a desktop environment with Python Tkinter installed.",
            file=sys.stderr,
        )
        sys.exit(1)
    EasyHapGUI().mainloop()


if __name__ == "__main__":
    main()
