# Supporting material

**[Paper title placeholder]** by [Author list placeholder]
***[Journal / special issue placeholder]***
[DOI/URL placeholdder]

## Overview

This repository contains the code, intermediate objects, and selected outputs used to support the analyses reported in the paper. The workflow combines R-based data preparation and post-processing with the `graph-tool` Python library to estimate a hierarchical stochastic blockmodel (SBM).

In broad terms, the repository covers four stages:

1. preparation of period-level FDI network data;
2. estimation of hierarchical block models with and without degree correction;
3. projection of hierarchical partitions back onto the original network;
4. post-processing, summary statistics, and figure generation.

## Repository contents

### Main scripts

- `01_Data.R`  
  Prepares the network data in R. It reads the raw input dataset, cleans country identifiers, constructs yearly and period-level matrices, computes the Balassa and logged-Balassa transformations, and exports descriptive outputs for later stages.

- `Paper.py` / `Paper_v2.py`  
  Python scripts for estimating hierarchical stochastic block models using `graph-tool`.

- `Input/export_hierarchy_partitions.py`  
  Utility script that converts pickled `graph-tool` nested blockmodel results into projected partition tables in CSV format.

- `02_Analysis.R`  
  Reads the projected partitions and prepared network objects, computes cluster-level summaries, processes World Bank indicators, and generates plots and analysis objects.

- `Paper.sh`  
  Example SLURM launcher for running the Python blockmodel stage inside an Apptainer container.

### Main folders

- `Input/`  
  Input objects and helper scripts used by the R and Python workflow.

- `Output/`  
  Saved model objects and exported intermediate results from the blockmodel estimation stage.

- `partitions/`  
  Projected hierarchical partitions exported from the `graph-tool` results.

### Additional archived artefacts

The branch also contains a number of saved analysis artefacts at top level, including `.RData`, `.RDS`, `.dput`, `.md`, and World Bank `.tab` files. These are preserved to document the empirical workflow and to simplify inspection of intermediate results.

## Data inputs

### FDI dataset

`01_Data.R` expects a raw file named:

```text
Code/Input/dataset.csv
```

This dataset is used to construct the bilateral FDI network. The script comments indicate that the intended extraction is based on inward direct investment liabilities, debt instruments, all entities, no aggregate countries, over the 2011-2023 period.

### World Bank indicators

`02_Analysis.R` expects World Bank indicator tables with filenames ending in `_WB.tab`. In the original workflow these are read from:

```text
./Analysis/
```

Examples included in this branch are:

- `GDPppp_WB.tab`
- `GNIpp_WB.tab`
- `pop_WB.tab`

## Important note on paths

The scripts were clearly written inside a larger working project rather than against the GitHub branch structure alone. In particular:

- `01_Data.R` reads from and writes to `./Code/Input/` and `./Code/Output/`
- `02_Analysis.R` reads from `./Code/Input/`, `./Code/Output/`, and `./Code/partitions/`
- `02_Analysis.R` also reads from and writes to `./Analysis/` and `./Figures/`
- `Paper.py` and `Paper_v2.py` work with `./Input/` and `./Output/`

So, if you want to rerun the scripts unmodified, the safest approach is to place this repository inside a project directory as the `Code/` subfolder and create sibling `Analysis/` and `Figures/` directories.

A practical layout is:

```text
project_root/
├── Code/          # clone or copy this repository here
├── Analysis/      # R analysis artefacts and *_WB.tab files
└── Figures/       # exported figures
```

If you prefer to keep the repository exactly as cloned, you will need to edit the hard-coded paths in the scripts.

## Software requirements

### R

A recent R installation is required for the two R scripts. The code uses packages from the following families:

- data manipulation and reshaping;
- plotting;
- country-code conversion;
- network and blockmodel utilities;
- table generation.

At minimum, inspect `01_Data.R` and `02_Analysis.R` before running them and install any missing packages manually.

### Python

The Python stage requires:

- Python 3
- `graph-tool`
- `numpy`

The helper script in `Input/export_hierarchy_partitions.py` also uses standard-library modules plus `numpy` and `graph-tool`.

### HPC / container execution

`Paper.sh` is a SLURM submission script and assumes an Apptainer image containing `graph-tool`. It is included as an example of the execution environment used by the author.

## Reproducing the workflow

## 1. Prepare the R inputs

Place the raw FDI dataset at:

```text
Code/Input/dataset.csv
```

Then run:

```bash
Rscript 01_Data.R
```

This stage generates the main network objects used downstream, including period-level matrices and logged-Balassa transformations.

## 2. Estimate hierarchical block models

Run one of the Python scripts:

```bash
python3 Paper_v2.py
```

or, on an HPC system configured like the original environment:

```bash
sbatch Paper.sh
```

The Python stage writes pickled blockmodel results to `Output/`.

## 3. Export projected partitions

Convert the pickled `graph-tool` results into CSV partitions. For example:

```bash
python3 Input/export_hierarchy_partitions.py \
  ./Output/res_HSBM-withDC-MCMC.pkl \
  ./Output/res_HSBM-woutDC-MCMC.pkl \
  --output-dir ./partitions \
  --node-prop name
```

This produces projected hierarchy-level assignments for each node.

## 4. Run the post-processing analysis

Ensure that:

- the projected partition CSV files are available where `02_Analysis.R` expects them;
- the World Bank `_WB.tab` files are available in `./Analysis/` if you are using the original paths;
- any required saved layout or auxiliary files are located consistently with the script paths.

Then run:

```bash
Rscript 02_Analysis.R
```

This stage produces:

- cluster summaries;
- World Bank summary tables by cluster and period;
- Sankey plots;
- mesoscopic network plots;
- mesoscopic matrix plots;
- saved `.RData` / `.RDS` analysis objects.

## Outputs

The repository includes or generates outputs such as:

- projected partition CSVs in `partitions/`;
- pickled `graph-tool` model objects in `Output/`;
- cluster assignment objects and ordering files;
- descriptive statistics and zero-value diagnostics;
- cluster profile markdown tables;
- mesoscopic graph objects and matrix objects;
- figure files exported to `Figures/`.

## Reproducibility remarks

This branch is best understood as an archival research-companion repository rather than a turnkey package. It preserves the main scripts and many intermediate objects, but some execution paths remain tied to the author's local project structure.

Accordingly, full reruns may require one or more of the following:

- recreating the original directory structure;
- moving archived files into the locations expected by the scripts;
- editing hard-coded paths;
- supplying the raw input dataset;
- installing the exact computational environment used for the `graph-tool` stage.

## Citation

> [Citation placeholdder].

## Licence

This repository is distributed under the Creative Common Attribution-NonCommercial-ShareAlike 4.0 International licence. See `LICENSE` for details.
