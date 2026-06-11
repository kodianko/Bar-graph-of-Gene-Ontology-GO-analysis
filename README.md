# Python scripts for proteomics data analysis, including Gene Ontology enrichment plots, volcano plots, heatmaps, and nanoparticle corona/lysosomal proteomics workflows.
Python code (.py) for generating bar graphs for proteomic data. Example processed proteomic data (.xlsx) are included.

## How to plot graphs

1. Download and extract zip file containing Data file (xlsx) and Python file (py):
 - NP_coronas_data.xlsx
 - Lysosomal_proteomics_data.xlsx
 - NP_coronas_proteins_groups_for_heatmap.xlsx
 - GO_BarPlots_NPCorona.py
 - VolcanoPlot_NP_coronas.py
 - VolcanoPlot_Lysosomal_proteomics.py
 - VolcanoPlot_Lysosomal_proteomics_groups.py
 - Heatmap_NP_coronas.py
 - Heatmap_Lysosomal_proteomics.py
 - Heatmap_Lysosomal_proteomics_Transp_Endo_Auto.py
 - Heatmap_Lysosomal_proteomics_Ferroptosis.py
 - Heatmap_Lysosomal_proteomics_Exo.py

2. Open Python file 
   We used PyCharm 2025.3.4

Requirements (listed in py files):
pandas
numpy
seaborn
matplotlib.pylab
matplotlib.colors

3. Modify the xlsx file path "df_path" to the location where you saved the downloaded xlsx file.


4. Change the "save_to_folder" setting where bar graphs will be stored:
 - replace 'X' with the name of the disk and its folders.

5. Run the code.

Cite the code: [![DOI](https://zenodo.org/badge/678244360.svg)](https://zenodo.org/badge/latestdoi/678244360)
