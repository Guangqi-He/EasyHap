from __future__ import annotations

import os
import threading
import tkinter as tk
from tkinter import colorchooser, filedialog, messagebox, ttk
from typing import Optional

from . import __version__
from .core import run_analysis


class EasyHapGUI(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title(f"EasyHap {__version__} — Research Edition")
        self.geometry("1000x820")
        self.minsize(850, 650)
        self._style()
        self._build()

    def _style(self):
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure("Title.TLabel", font=("Arial", 18, "bold"))
        style.configure("Sub.TLabel", font=("Arial", 10))
        style.configure("Req.TLabel", font=("Arial", 10, "bold"))
        style.configure("Run.TButton", font=("Arial", 11, "bold"), padding=(14, 8))
        style.configure("TLabelframe.Label", font=("Arial", 10, "bold"))

    def _file_row(self, parent, row, label, var, required=False, directory=False):
        ttk.Label(parent, text=label + ("  *" if required else "  (optional)"), style="Req.TLabel" if required else "TLabel").grid(row=row,column=0,sticky="w",padx=8,pady=5)
        ttk.Entry(parent,textvariable=var).grid(row=row,column=1,sticky="ew",padx=8,pady=5)
        def choose():
            path=filedialog.askdirectory() if directory else filedialog.askopenfilename()
            if path: var.set(path)
        ttk.Button(parent,text="Browse",command=choose).grid(row=row,column=2,padx=8,pady=5)

    def _build(self):
        root=ttk.Frame(self,padding=14); root.pack(fill="both",expand=True)
        ttk.Label(root,text="EasyHap",style="Title.TLabel").pack(anchor="w")
        ttk.Label(root,text="Ploidy-aware haplotype analysis, population comparison, trait association, LD and publication-ready visualization",style="Sub.TLabel").pack(anchor="w",pady=(0,10))
        ttk.Label(root,text="* Required. Provide a VCF plus either a single Region or a Region file. Other inputs are optional.").pack(anchor="w",pady=(0,8))

        nb=ttk.Notebook(root); nb.pack(fill="both",expand=True)
        inp=ttk.Frame(nb,padding=12); ana=ttk.Frame(nb,padding=12); vis=ttk.Frame(nb,padding=12)
        nb.add(inp,text="Inputs"); nb.add(ana,text="Analysis"); nb.add(vis,text="Visualization")
        inp.columnconfigure(1,weight=1); ana.columnconfigure(1,weight=1); vis.columnconfigure(1,weight=1)

        self.vcf=tk.StringVar(); self.group=tk.StringVar(); self.traits=tk.StringVar(); self.gff=tk.StringVar(); self.outdir=tk.StringVar(value="EasyHap_results")
        self.region=tk.StringVar(); self.region_file=tk.StringVar()
        self._file_row(inp,0,"Phased VCF / VCF.GZ / BCF",self.vcf,True)
        ttk.Label(inp,text="Region  *",style="Req.TLabel").grid(row=1,column=0,sticky="w",padx=8,pady=5)
        ttk.Entry(inp,textvariable=self.region).grid(row=1,column=1,sticky="ew",padx=8,pady=5)
        ttk.Label(inp,text="e.g. Chr10:1-500").grid(row=1,column=2,sticky="w",padx=8)
        self._file_row(inp,2,"Region file (batch)",self.region_file,False)
        ttk.Label(inp,text="Use either Region or Region file.").grid(row=3,column=1,sticky="w",padx=8,pady=(0,8))
        self._file_row(inp,4,"Sample-group file (optional lat/lon in columns 3-4)",self.group,False)
        self._file_row(inp,5,"Trait table",self.traits,False)
        self._file_row(inp,6,"Gene annotation (GFF3/GTF)",self.gff,False)
        self._file_row(inp,7,"Output directory",self.outdir,False,True)

        self.mode=tk.StringVar(value="inbred"); self.hetero=tk.StringVar(value="slash"); self.min_variants=tk.StringVar(value="2")
        self.cluster=tk.StringVar(value="0.15"); self.fisher_groups=tk.StringVar(); self.fisher_alpha=tk.StringVar(); self.fisher_adjust=tk.StringVar(value="none")
        self.ld=tk.BooleanVar(value=True); self.processed=tk.BooleanVar(value=True); self.network=tk.BooleanVar(value=True)
        self.gene_feature=tk.StringVar(value="all")
        ttk.Label(ana,text="Mode").grid(row=0,column=0,sticky="w",padx=8,pady=6); ttk.Combobox(ana,textvariable=self.mode,values=["inbred","hybrid"],state="readonly",width=14).grid(row=0,column=1,sticky="w")
        ttk.Label(ana,text="Heterozygous-site encoding (inbred mode)").grid(row=1,column=0,sticky="w",padx=8,pady=6); ttk.Combobox(ana,textvariable=self.hetero,values=["slash","iupac","missing"],state="readonly",width=14).grid(row=1,column=1,sticky="w")
        ttk.Label(ana,text="Minimum variants per region").grid(row=2,column=0,sticky="w",padx=8,pady=6); ttk.Entry(ana,textvariable=self.min_variants,width=10).grid(row=2,column=1,sticky="w")
        ttk.Label(ana,text="Haplotype cluster distance threshold").grid(row=3,column=0,sticky="w",padx=8,pady=6); ttk.Entry(ana,textvariable=self.cluster,width=10).grid(row=3,column=1,sticky="w")
        ttk.Separator(ana).grid(row=4,column=0,columnspan=3,sticky="ew",pady=10)
        ttk.Label(ana,text="Fisher filtering groups").grid(row=5,column=0,sticky="w",padx=8,pady=6); ttk.Entry(ana,textvariable=self.fisher_groups,width=28).grid(row=5,column=1,sticky="w"); ttk.Label(ana,text="e.g. Cultivar,Landrace").grid(row=5,column=2,sticky="w")
        ttk.Label(ana,text="Fisher significance threshold").grid(row=6,column=0,sticky="w",padx=8,pady=6); ttk.Entry(ana,textvariable=self.fisher_alpha,width=10).grid(row=6,column=1,sticky="w")
        ttk.Label(ana,text="Multiple-testing adjustment").grid(row=7,column=0,sticky="w",padx=8,pady=6); ttk.Combobox(ana,textvariable=self.fisher_adjust,values=["none","bh"],state="readonly",width=10).grid(row=7,column=1,sticky="w")
        ttk.Checkbutton(ana,text="Calculate LD r² and LD heatmap",variable=self.ld).grid(row=8,column=0,sticky="w",padx=8,pady=6)
        ttk.Checkbutton(ana,text="Write processed allele/genotype token tables",variable=self.processed).grid(row=9,column=0,columnspan=2,sticky="w",padx=8,pady=6)
        ttk.Label(ana,text="Gene feature used for variant analysis").grid(row=10,column=0,sticky="w",padx=8,pady=6)
        ttk.Combobox(ana,textvariable=self.gene_feature,values=["all","exon","intron","cds","utr"],state="readonly",width=14).grid(row=10,column=1,sticky="w")
        ttk.Label(ana,text="exon/intron/CDS/UTR requires a GFF3/GTF file").grid(row=10,column=2,sticky="w")

        self.plot=tk.BooleanVar(value=True); self.plot_format=tk.StringVar(value="pdf"); self.level=tk.StringVar(value="hap"); self.plot_min=tk.StringVar(value="1"); self.trait_cols=tk.StringVar()
        self.ref_color=tk.StringVar(value="#70AD47"); self.alt_color=tk.StringVar(value="#4472C4"); self.missing_color=tk.StringVar(value="#D9D9D9")
        self.hap_palette=tk.StringVar(value="#4E79A7,#F28E2B,#E15759,#76B7B2,#59A14F,#EDC948,#B07AA1,#FF9DA7,#9C755F,#BAB0AC")
        self.ld_cmap=tk.StringVar(value="viridis")
        self.cell_text=tk.StringVar(value="auto")
        self.map_style=tk.StringVar(value="auto")
        ttk.Checkbutton(vis,text="Generate publication figures",variable=self.plot).grid(row=0,column=0,columnspan=2,sticky="w",padx=8,pady=6)
        ttk.Label(vis,text="Figure formats").grid(row=1,column=0,sticky="w",padx=8,pady=6); ttk.Entry(vis,textvariable=self.plot_format,width=18).grid(row=1,column=1,sticky="w")
        ttk.Label(vis,text="Displayed haplotype level").grid(row=2,column=0,sticky="w",padx=8,pady=6); ttk.Combobox(vis,textvariable=self.level,values=["hap","cluster"],state="readonly",width=12).grid(row=2,column=1,sticky="w")
        ttk.Label(vis,text="Minimum displayed haplotype/cluster count").grid(row=3,column=0,sticky="w",padx=8,pady=6); ttk.Entry(vis,textvariable=self.plot_min,width=10).grid(row=3,column=1,sticky="w")
        ttk.Label(vis,text="Trait columns to plot").grid(row=4,column=0,sticky="w",padx=8,pady=6); ttk.Entry(vis,textvariable=self.trait_cols).grid(row=4,column=1,sticky="ew"); ttk.Label(vis,text="comma-separated; blank = all").grid(row=4,column=2,sticky="w")
        ttk.Separator(vis).grid(row=5,column=0,columnspan=3,sticky="ew",pady=10)
        for r,(lab,var) in enumerate([("Haplotype heatmap: REF color",self.ref_color),("Haplotype heatmap: ALT color",self.alt_color),("Haplotype heatmap: missing color",self.missing_color)],6):
            ttk.Label(vis,text=lab).grid(row=r,column=0,sticky="w",padx=8,pady=6); ttk.Entry(vis,textvariable=var,width=14).grid(row=r,column=1,sticky="w")
            ttk.Button(vis,text="Palette…",command=lambda v=var:self._pick_color(v)).grid(row=r,column=2,sticky="w",padx=8)
        ttk.Label(vis,text="Haplotype color palette").grid(row=9,column=0,sticky="w",padx=8,pady=6)
        ttk.Entry(vis,textvariable=self.hap_palette).grid(row=9,column=1,sticky="ew")
        ttk.Label(vis,text="comma-separated; shared by pie, stacked bar and boxplot").grid(row=9,column=2,sticky="w")
        ttk.Label(vis,text="LD heatmap palette").grid(row=10,column=0,sticky="w",padx=8,pady=6)
        ttk.Entry(vis,textvariable=self.ld_cmap,width=18).grid(row=10,column=1,sticky="w")
        ttk.Label(vis,text="Matplotlib colormap, e.g. viridis, magma, coolwarm").grid(row=10,column=2,sticky="w")
        ttk.Label(vis,text="Allele text in heatmap cells").grid(row=11,column=0,sticky="w",padx=8,pady=6)
        ttk.Combobox(vis,textvariable=self.cell_text,values=["auto","always","never"],state="readonly",width=12).grid(row=11,column=1,sticky="w")
        ttk.Label(vis,text="auto = dynamic font; always = force labels; never = hide labels").grid(row=11,column=2,sticky="w")
        ttk.Label(vis,text="Geographic haplotype map").grid(row=12,column=0,sticky="w",padx=8,pady=6)
        ttk.Combobox(vis,textvariable=self.map_style,values=["auto","pie","bar","both","none"],state="readonly",width=12).grid(row=12,column=1,sticky="w")
        ttk.Label(vis,text="auto = pie when sample-group latitude/longitude are available").grid(row=12,column=2,sticky="w")
        ttk.Checkbutton(vis,text="Generate minimum-spanning haplotype network",variable=self.network).grid(row=13,column=0,columnspan=2,sticky="w",padx=8,pady=6)

        bottom=ttk.Frame(root); bottom.pack(fill="x",pady=(10,0))
        ttk.Button(bottom,text="Run EasyHap",style="Run.TButton",command=self._run).pack(side="left")
        self.status=tk.StringVar(value="Ready"); ttk.Label(bottom,textvariable=self.status).pack(side="left",padx=14)
        self.log=tk.Text(root,height=10,font=("Consolas",9)); self.log.pack(fill="x",pady=(8,0))

    def _pick_color(self,var):
        color=colorchooser.askcolor(color=var.get(),title="Choose allele color")[1]
        if color: var.set(color)

    def _log(self,text):
        self.log.insert("end",text+"\n"); self.log.see("end")

    def _run(self):
        if not self.vcf.get().strip(): messagebox.showerror("Missing required input","VCF is required."); return
        if bool(self.region.get().strip()) == bool(self.region_file.get().strip()): messagebox.showerror("Region input","Provide exactly one of Region or Region file."); return
        try:
            minv=int(self.min_variants.get()); plotmin=int(self.plot_min.get()); cluster=float(self.cluster.get())
            if minv<1 or plotmin<1: raise ValueError
        except ValueError: messagebox.showerror("Invalid numeric option","Minimum variants/count must be positive integers and cluster threshold numeric."); return
        fg1=fg2=None; alpha=None
        if self.fisher_groups.get().strip():
            parts=[x.strip() for x in self.fisher_groups.get().split(",") if x.strip()]
            if len(parts)!=2 or not self.group.get().strip(): messagebox.showerror("Fisher filtering","Fisher filtering requires a group file and exactly two group names."); return
            fg1,fg2=parts
        if self.fisher_alpha.get().strip():
            try: alpha=float(self.fisher_alpha.get())
            except ValueError: messagebox.showerror("Fisher significance threshold","Fisher alpha must be numeric."); return
        if self.gene_feature.get() != "all" and not self.gff.get().strip():
            messagebox.showerror("Gene feature filter", f"Gene feature '{self.gene_feature.get()}' requires a GFF3/GTF file."); return
        out=self.outdir.get().strip() or "EasyHap_results"
        if self.plot.get():
            from .plotting import set_cell_text_mode
            set_cell_text_mode(self.cell_text.get())
        kwargs=dict(vcf_path=self.vcf.get().strip(),group_file=self.group.get().strip() or None,outdir=out,region=self.region.get().strip() or None,region_file=self.region_file.get().strip() or None,
                    mode=self.mode.get(),hetero_policy=self.hetero.get(),trait_file=self.traits.get().strip() or None,gff_file=self.gff.get().strip() or None,gene_feature=self.gene_feature.get(),
                    fisher_group1=fg1,fisher_group2=fg2,fisher_alpha=alpha,fisher_adjust=self.fisher_adjust.get(),cluster_threshold=cluster,min_variants=minv,
                    write_processed=self.processed.get(),make_plots=self.plot.get(),plot_formats=[x.strip() for x in self.plot_format.get().split(",") if x.strip()],
                    traits_to_plot=[x.strip() for x in self.trait_cols.get().split(",") if x.strip()] or None,plot_hap_level=self.level.get(),plot_min_count=plotmin,
                    ref_color=self.ref_color.get(),alt_color=self.alt_color.get(),missing_color=self.missing_color.get(),
                    hap_palette=[x.strip() for x in self.hap_palette.get().split(",") if x.strip()] or None,ld_cmap=self.ld_cmap.get().strip() or "viridis",make_ld=self.ld.get(),
                    make_network=self.network.get(),map_style=self.map_style.get())
        self.status.set("Running…"); self._log(f"Running EasyHap {__version__}…")
        def worker():
            try:
                res=run_analysis(**kwargs)
                self.after(0,lambda:self.status.set(f"Finished: {len(res)} region(s)")); self.after(0,lambda:self._log(f"Finished {len(res)} region(s). Log: {os.path.join(out,'EasyHap.log')}"))
            except Exception as exc:
                self.after(0,lambda:self.status.set("Error")); self.after(0,lambda:self._log(f"ERROR: {exc}")); self.after(0,lambda:messagebox.showerror("EasyHap error",str(exc)))
        threading.Thread(target=worker,daemon=True).start()


def main():
    EasyHapGUI().mainloop()

if __name__ == "__main__": main()
