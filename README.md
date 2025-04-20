# PURIFI
SNP-Based QA/QC Pipeline Shiny App
This repository contains an **R Shiny** application for **genetic purity** analysis and **F1 pedigree verification** using SNP marker data, tailored for dryland crop breeding programs.

# SNP-Based QA/QC Pipeline Shiny App

---

## 📁 Repository Structure

```text
├── app.R                 # Main Shiny application (ui + server)
├── functions/            # Modular R scripts with core functions
│   ├── read_data.R
│   ├── grid_preparation.R
│   ├── snp_summary.R
│   ├── parent_consensus.R
│   ├── purity_calculation.R
│   ├── f1_verification.R
│   ├── hapmap_export.R
├── www/                  # Static assets (example data, images, CSS)
│   ├── meta.txt          # Example metadata
│   ├── geno.csv          # Example genotype file
│   ├── logoW.png
│   └── img*.png
├── README.md             # Project documentation (this file)
├── LICENSE
└── .gitignore
```

---

## 🚀 Getting Started

### Prerequisites

- R (>= 4.0)
- RStudio (optional but recommended)
- Internet connection to install packages

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/YourOrg/snp-qaqc-pipeline.git
   cd snp-qaqc-pipeline
   ```
2. **Install required R packages**
   ```r
   install.packages(c(
     "shiny", "shinythemes", "shinyjs", 
     "readxl", "tidyverse", "DT", "gt", 
     "plotly", "webshot", "shinycssloaders"
   ))
   ```
3. **Launch the Shiny App**
   ```r
   library(shiny)
   runApp("app.R")
   ```

---

## 🖥️ Application Workflow

1. **Data Upload**

   - Upload sample metadata (.txt) and SNP genotype data (.csv) via the UI.
   - Optionally load example datasets from `www/`.

2. **Marker Visualization & Summary**

   - Interactive SNP cluster plots (X vs Y) to assess call quality.
   - Summary table with missing call rates, allele frequencies, PIC, etc.

3. **Parent Consensus & QC**

   - Build consensus genotypes for parent lines.
   - Calculate purity scores and visualize parent sample consistency.
   - Adjustable purity threshold slider.

4. **F1 Verification**

   - Compare F1 genotypes against parent consensus.
   - Calculate heterozygosity percentage and classify F1 crosses.
   - Customizable heterozygosity threshold slider.

5. **HAPMAP Export**

   - Clean genotype data exported in HAPMAP format for downstream tools.

---

## 🔍 Key Functions (in `functions/`)

| File                   | Function(s)                                        |
| ---------------------- | -------------------------------------------------- |
| `read_data.R`          | `read_datasets()`                                  |
| `grid_preparation.R`   | `create_grid_file()`                               |
| `snp_summary.R`        | `calculate_snp_info()`, `compute_sample_summary()` |
| `parent_consensus.R`   | `create_parents_consensus()`                       |
| `purity_calculation.R` | `parental_purity()`, `create_purity_report()`      |
| `f1_verification.R`    | `find_homo_poly_snps()`, `create_f1_report()`      |
| `hapmap_export.R`      | `create_hapmap()`                                  |

---

## 🎨 Customization & Configuration

- **Thresholds** in the UI:
  - Parent Purity (default: **80%**)
  - F1 Heterozygosity (default: **60%**)
- **Styling**: Modify CSS in `app.R` head tags for color schemes and table style.
- **Example Data**: Located in `www/meta.txt` and `www/geno.csv`.

---

## 🚧 Deployment

- **Local**: `runApp("app.R")` in RStudio or R console.
- **Shiny Server**: Place this directory under `/srv/shiny-server/snp-qaqc-pipeline`.
- **shinyapps.io**:
  ```bash
  rsconnect::deployApp()
  ```

## Contributors

- **Dr. Abhishek Rathore** — Principal Scientist, Breeding Data & Informatics
- **Mr. Peter Kimathi** — Bioinformatics & Software Developer
- **Dr. Roma Rani Das** — Biometrician

---

## 📄 License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

