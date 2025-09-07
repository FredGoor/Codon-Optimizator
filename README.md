# Codon Optimizer (MATLAB)

This repository contains a MATLAB pipeline for **codon usage analysis and optimization** of coding DNA sequences (CDS).  
The pipeline allows you to compare codon usage of an input gene against reference codon profiles and to generate a recoded version of the gene optimized toward a specified codon usage bias (e.g. ribosomal genes, horizontally transferred genes, or custom codon tables).  

---

## ✨ Features

- **Codon usage analysis**
  - Extract GC content, absolute codon usage, and relative synonymous codon usage (RSU) from any input CDS (FASTA/TXT).
  - Export results into Excel spreadsheets.

- **Codon usage visualization**
  - Bar plots of codon usage.
  - Scatter plots with regression and R².
  - Ratio plots comparing two codon usage profiles.

- **Codon recoding**
  - Provide a codon usage table (reference bias) to guide recoding.
  - Generate a synonymous DNA sequence preserving the protein sequence but matching the target codon usage.
  - Verify amino acid integrity (original vs optimized sequence).

- **Output**
  - Multi-sheet Excel file including:
    - Original gene codon usage
    - Protein sequence
    - Reference codon usage table
    - Optimized DNA sequence
    - Optimized codon usage
    - Comparison plots

---

## 📂 Repository structure

- `CodonOptimizer_main.m` – main pipeline script.  
- `Fun_Fred_CodonAnalysis.m` – analyze codon usage and GC content from a CDS.  
- `Fun_Fred_CodonAnalysis_bis.m` – analyze codon usage from a pre-computed codon table.  
- `Fun_Fred_Codontable_comparison.m` – compare two codon usage tables (bar + scatter plots).  
- `Fun_Fred_scatter_reg.m` – scatter plot with regression line and R².  
- `Fun_Fred_Reg_Linear.m` – linear regression helper.  
- `Fun_Fred_Codongraphs.m` – bar graph visualization of codon usage.
- 'Example folder' containing mNeonGreen cds and resulting optimized sequence mNeonGreen_Optiseq.txt

---

## 📖 Example

The `example/` folder provides test files:  

- `mNeonGreen.txt` – input CDS (wild-type sequence).  
- `ribosomal_codon_table.txt` – example codon usage bias (ribosomal genes).  
- `mNeonGreen_optiseq.txt` – optimized version of mNeonGreen generated with the pipeline.  


## 🚀 Usage

1. Open `CodonOptimizer_main.m` in MATLAB.  
2. Run the script and select:
   - The CDS file (FASTA or TXT, single coding sequence).
   - The target codon usage table (TXT, 2 columns: codon / frequency or AA-codon / frequency).  
3. Results will be saved as an Excel file in the same folder as the input sequence.  

---

## 📖 Example

- Input CDS: `example_gene.fasta`  
- Target codon usage: `ribosomal_codon_table.txt`  
- Output: `example_gene_datatable.xlsx` with codon usage, optimized DNA sequence, protein sequence, and plots.  

---

## 🔧 Requirements

- MATLAB (tested ≥ R2021b).  
- Bioinformatics Toolbox (for `fastaread`, `nt2aa`, `codoncount`).  

---

## 📜 License

MIT License. Feel free to use, modify, and share.  
