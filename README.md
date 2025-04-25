# PURIFI
SNP-Based QA/QC Pipeline Shiny App
This repository contains an **R Shiny** application for **genetic purity** analysis and **F1 pedigree verification** using SNP marker data, tailored for Dryland Crop Programs (DCP) of CIMMYT.

# SNP-Based QA/QC Pipeline Shiny App

---
# PURIFI: SNP-Based QA/QC Pipeline Shiny Application 

This repository contains a single-file **R Shiny** application (`app.R`) called **PURIFI**, designed for **genetic purity** analysis and **F1 pedigree verification** using SNP marker data, tailored for dryland crop breeding programs.

---

## 📁 Repository Structure

```text
├── app.R             # Main Shiny application (UI and server logic)
├── www/              # Static assets (example data, images, CSS)
│   ├── meta.txt      # Example sample metadata
│   ├── geno.csv      # Example genotype data
│   ├── logoW.png     # Logo image
│   └── img1.png, img2.png, img3.png  # Panel icons
├── README.md         # Project documentation (this file)
├── LICENSE           # MIT License
└── .gitignore        # Files and folders to ignore in Git
```

---

## 🚀 Getting Started

### Prerequisites
- R (>= 4.0)
- RStudio (optional but recommended)

### Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/reach2abhi/PURIFI.git
   cd PURIFI
   ```
2. **Install required R packages**:
   ```r
   install.packages(c(
     "shiny", "shinythemes", "shinyjs", 
     "readxl", "tidyverse", "DT", "gt", 
     "plotly", "webshot", "shinycssloaders"
   ))
   ```
3. **Launch the Shiny App**:
   ```r
   library(shiny)
   runApp("app.R")
   ```

---

## 🖥️ Application Workflow

1. **Data Upload**
   - Upload sample metadata (.txt) and SNP genotype data (.csv) through the UI.
   - Load provided example files from the `www/` folder.

2. **Marker Visualization & Summary**
   - Interactive SNP cluster plots (X vs Y) to evaluate call quality.
   - Summary table of each marker’s missing data rate, allele frequencies, heterozygosity, and PIC.

3. **Parent Consensus & QC**
   - Generate consensus genotypes across parent replicates.
   - Compute parent purity scores with a customizable threshold slider.
   - Visual heatmap of parent sample consistency.

4. **F1 Verification**
   - Compare F1 hybrid genotypes to parent consensus.
   - Calculate percentage heterozygosity and classify crosses as:
     - **Successful F1** (100% heterozygous)
     - **Possible F1** (≥ threshold and <100%)
     - **Parent Quality Failed** (purity score below threshold)
     - **Failed F1** (heterozygosity below threshold)
   - Customizable heterozygosity threshold slider.

5. **HAPMAP Export**
   - Export cleaned genotype data into a HAPMAP-format file (`hapmap-<date>.hmp.txt`) for downstream analysis.

---

## 🔍 Core Functions (defined in `app.R`)

| Function Name            | Purpose                                                     |
| ------------------------ | ----------------------------------------------------------- |
| `read_datasets()`        | Reads metadata and raw genotype files, parses SNP calls.   |
| `create_grid_file()`     | Reshapes raw calls into a Sample×SNP matrix.               |
| `calculate_snp_info()`   | Computes marker-level statistics (missing %, PIC, etc.).    |
| `create_parents_consensus()` | Builds consensus genotypes for parent lines.           |
| `parental_purity()`      | Calculates purity comparisons between samples and consensus.|
| `create_f1_report()`     | Evaluates F1 hybrids for expected heterozygosity.          |
| `create_hapmap()`        | Formats final data into HAPMAP standard output.            |

---

## 🎨 Customization & Settings

- **Threshold sliders** in the UI:
  - Parent Purity (default **80%**)
  - F1 Heterozygosity (default **60%**)
- **Styling/CSS** tweaks in the `<head>` section of `app.R`.
- **Example data** in `www/` can be replaced with project-specific files.

---

## 🚧 Deployment

- **Locally**: `shiny::runApp("app.R")`.
- **Shiny Server**: Copy this folder to `/srv/shiny-server/PURIFI/`.
- **shinyapps.io**: Use `rsconnect::deployApp()` from within the project directory.

---

## 👥 Contributors

- **Dr. Abhishek Rathore** — Principal Scientist, Breeding Data & Informatics
- **Mr. Peter Kimathi** — Bioinformatics & Software Developer
- **Dr. Roma Rani Das** — Biometrician

---

## 📄 License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.



