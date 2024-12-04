# Digital Mental Health Interventions for Depression: A Multiverse Meta-Analysis

This repository contains the data, code, and resources used for the publication **"Digital mental health interventions for the treatment of depression: A multiverse meta-analysis"**, available on [PubMed](https://pubmed.ncbi.nlm.nih.gov/39419189/).

## Repository Structure

### Files and Folders
- **`data/`**  
  Contains the datasets used for the analysis. These datasets include study-level data on digital interventions for depression, formatted for use in the R scripts.

- **`functions/`**  
  Includes helper functions used in the analysis and manuscript preparation. These functions handle tasks such as data processing, statistical modeling, and figure generation.

- **`digDep_1_multiverse.Rmd`**  
  This R Markdown file sets up the multiverse meta-analysis. It runs multiple analytic models to examine the impact of analytic choices on the results.

- **`digDep_2_results.Rmd`**  
  This R Markdown file produces the main results of the study. It summarizes the findings from the multiverse meta-analysis, including tables and descriptive statistics.

- **`digDep_3_manuscript.Rmd`**  
  This R Markdown file generates the figures and tables used in the manuscript. It includes code for creating visualizations that illustrate the results of the multiverse analysis.

- **`digitalInterventionsDepression.Rproj`**  
  The RStudio project file for this repository. Open this file in RStudio to access the project environment.

- **`.gitignore`**  
  Specifies files and folders that are ignored by Git version control (e.g., temporary files or outputs).

- **`README.md`**  
  This file. Provides an overview of the repository structure and its purpose.

## How to Use
1. Download this repository as a ZIP file by clicking the green **"Code"** button at the top of the GitHub page and selecting **"Download ZIP"**. or clone this repository to your local machine:
```bash
git clone https://github.com/cyplessen/digital-interventions-depression
```
2. Extract the contents of the ZIP file to a folder on your computer.
3. Open `digitalInterventionsDepression.Rproj` in RStudio to load the project environment.
4. Run the R Markdown files in sequence to replicate the analysis and results:
   - Start with `digDep_1_multiverse.Rmd` to execute the multiverse analysis.
   - Use `digDep_2_results.Rmd` to generate the main results.
   - Run `digDep_3_manuscript.Rmd` to create figures and tables for the manuscript.

## Citation
If you use the code or data in this repository, please cite the publication:

> **Digital mental health interventions for the treatment of depression: A multiverse meta-analysis**  
> Constantin Yves Plessen*, Olga Maria Panagiotopoulou*, Lingyao Tong, Pim Cuijpers, Eirini Karyotaki
> * Shared first authors.
> Published in  J Affect Disord.  
> DOI: [https:://10.1016/j.jad.2024.10.018.](https://doi.org/10.1016/j.jad.2024.10.018)

## Contact
For questions or feedback about this repository, please contact me.
