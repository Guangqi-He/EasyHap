from __future__ import annotations

import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Optional

from .core import run_analysis


class EasyHapGUI(tk.Tk):
    """Minimal Windows-friendly GUI wrapper around the same Python core."""

    def __init__(self) -> None:
        super().__init__()
        self.title("EasyHap 1.0")
        self.geometry("760x520")
        self._build()

    def _row_file(self, parent, row: int, label: str, var: tk.StringVar, save_dir: bool = False) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=8, pady=4)
        ttk.Entry(parent, textvariable=var, width=72).grid(row=row, column=1, sticky="we", padx=8, pady=4)
        def choose():
            if save_dir:
                path = filedialog.askdirectory()
            else:
                path = filedialog.askopenfilename()
            if path:
                var.set(path)
        ttk.Button(parent, text="Browse", command=choose).grid(row=row, column=2, sticky="e", padx=8, pady=4)

    def _build(self) -> None:
        frm = ttk.Frame(self)
        frm.pack(fill="both", expand=True, padx=10, pady=10)
        frm.columnconfigure(1, weight=1)

        self.vcf = tk.StringVar()
        self.group = tk.StringVar()
        self.traits = tk.StringVar()
        self.gff = tk.StringVar()
        self.outdir = tk.StringVar()
        self.region = tk.StringVar(value="ChrDemo:900-2900")
        self.region_file = tk.StringVar()
        self.mode = tk.StringVar(value="inbred")
        self.hetero_policy = tk.StringVar(value="slash")
        self.fisher_groups = tk.StringVar()
        self.fisher_alpha = tk.StringVar()
        self.fisher_adjust = tk.StringVar(value="none")
        self.plot = tk.BooleanVar(value=True)
        self.plot_format = tk.StringVar(value="pdf")
        self.plot_hap_level = tk.StringVar(value="hap")
        self.plot_min_count = tk.StringVar(value="1")

        self._row_file(frm, 0, "VCF", self.vcf)
        self._row_file(frm, 1, "Sample group", self.group)
        self._row_file(frm, 2, "Trait table", self.traits)
        self._row_file(frm, 3, "GFF/GTF", self.gff)
        self._row_file(frm, 4, "Output directory", self.outdir, save_dir=True)

        ttk.Label(frm, text="Region").grid(row=5, column=0, sticky="w", padx=8, pady=4)
        ttk.Entry(frm, textvariable=self.region).grid(row=5, column=1, sticky="we", padx=8, pady=4)
        ttk.Label(frm, text="or region file").grid(row=6, column=0, sticky="w", padx=8, pady=4)
        ttk.Entry(frm, textvariable=self.region_file).grid(row=6, column=1, sticky="we", padx=8, pady=4)
        ttk.Button(frm, text="Browse", command=lambda: self._choose_file(self.region_file)).grid(row=6, column=2, padx=8, pady=4)

        opt = ttk.Frame(frm)
        opt.grid(row=7, column=0, columnspan=3, sticky="we", pady=8)
        ttk.Label(opt, text="Mode").grid(row=0, column=0, padx=8)
        ttk.Combobox(opt, textvariable=self.mode, values=["inbred", "hybrid"], width=12, state="readonly").grid(row=0, column=1)
        ttk.Label(opt, text="Hetero policy").grid(row=0, column=2, padx=8)
        ttk.Combobox(opt, textvariable=self.hetero_policy, values=["slash", "iupac", "missing"], width=12, state="readonly").grid(row=0, column=3)
        ttk.Checkbutton(opt, text="Plot", variable=self.plot).grid(row=0, column=4, padx=8)
        ttk.Label(opt, text="Format").grid(row=0, column=5, padx=8)
        ttk.Entry(opt, textvariable=self.plot_format, width=10).grid(row=0, column=6)
        ttk.Label(opt, text="Plot level").grid(row=1, column=0, padx=8, pady=4)
        ttk.Combobox(opt, textvariable=self.plot_hap_level, values=["hap", "cluster"], width=12, state="readonly").grid(row=1, column=1, pady=4)
        ttk.Label(opt, text="Plot min count").grid(row=1, column=2, padx=8, pady=4)
        ttk.Entry(opt, textvariable=self.plot_min_count, width=8).grid(row=1, column=3, pady=4)

        fisher = ttk.Frame(frm)
        fisher.grid(row=8, column=0, columnspan=3, sticky="we", pady=4)
        ttk.Label(fisher, text="Fisher groups, e.g. Cultivar,Landrace").grid(row=0, column=0, padx=8)
        ttk.Entry(fisher, textvariable=self.fisher_groups, width=24).grid(row=0, column=1)
        ttk.Label(fisher, text="alpha").grid(row=0, column=2, padx=8)
        ttk.Entry(fisher, textvariable=self.fisher_alpha, width=8).grid(row=0, column=3)
        ttk.Label(fisher, text="adjust").grid(row=0, column=4, padx=8)
        ttk.Combobox(fisher, textvariable=self.fisher_adjust, values=["none", "bh"], width=8, state="readonly").grid(row=0, column=5)

        ttk.Button(frm, text="Run", command=self._run).grid(row=9, column=0, sticky="w", padx=8, pady=8)
        self.log = tk.Text(frm, height=12)
        self.log.grid(row=10, column=0, columnspan=3, sticky="nsew", padx=8, pady=8)
        frm.rowconfigure(10, weight=1)

    def _choose_file(self, var: tk.StringVar) -> None:
        path = filedialog.askopenfilename()
        if path:
            var.set(path)

    def _log(self, text: str) -> None:
        self.log.insert("end", text + "\n")
        self.log.see("end")

    def _run(self) -> None:
        if not self.vcf.get() or not self.group.get() or not self.outdir.get():
            messagebox.showerror("Missing input", "VCF, sample group, and output directory are required.")
            return
        fisher_group1 = fisher_group2 = None
        if self.fisher_groups.get().strip():
            parts = [x.strip() for x in self.fisher_groups.get().split(",") if x.strip()]
            if len(parts) != 2:
                messagebox.showerror("Invalid Fisher groups", "Use exactly two groups, e.g. Cultivar,Landrace")
                return
            fisher_group1, fisher_group2 = parts
        alpha: Optional[float] = None
        if self.fisher_alpha.get().strip():
            try:
                alpha = float(self.fisher_alpha.get())
            except ValueError:
                messagebox.showerror("Invalid alpha", "Fisher alpha should be a number.")
                return

        try:
            plot_min_count = int(self.plot_min_count.get().strip() or "1")
            if plot_min_count < 1:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid plot min count", "Plot min count should be a positive integer.")
            return

        kwargs = dict(
            vcf_path=self.vcf.get(),
            group_file=self.group.get(),
            outdir=self.outdir.get(),
            region=self.region.get().strip() or None,
            region_file=self.region_file.get().strip() or None,
            mode=self.mode.get(),
            hetero_policy=self.hetero_policy.get(),
            trait_file=self.traits.get().strip() or None,
            fisher_group1=fisher_group1,
            fisher_group2=fisher_group2,
            fisher_alpha=alpha,
            fisher_adjust=self.fisher_adjust.get(),
            make_plots=self.plot.get(),
            gff_file=self.gff.get().strip() or None,
            plot_formats=[x.strip() for x in self.plot_format.get().split(",") if x.strip()],
            plot_hap_level=self.plot_hap_level.get(),
            plot_min_count=plot_min_count,
        )
        self._log("Running...")
        def worker():
            try:
                results = run_analysis(**kwargs)
                self.after(0, lambda: self._log(f"Finished {len(results)} region(s)."))
                for r in results:
                    self.after(0, lambda r=r: self._log(f"{r.region.vcf_label}: {r.output_prefix}"))
            except Exception as exc:
                self.after(0, lambda: messagebox.showerror("EasyHap error", str(exc)))
                self.after(0, lambda: self._log(f"ERROR: {exc}"))
        threading.Thread(target=worker, daemon=True).start()


def main() -> None:
    app = EasyHapGUI()
    app.mainloop()


if __name__ == "__main__":
    main()
