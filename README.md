# Connectional_Variability_Axis
Data and codes for our paper **"Connectional Axis of Individual Functional Variability: Patterns, structural correlates, and implications for development and cognition**.

*Inter-individual fc variability* matrices for the **HCP-D** and **HCP-YA** (typically referred as **'HCP'**) datasets, and the connectome matrices of *functional connectivity* and *individual variability in structural communicability* are obtained using `Schaefer-400` (7 Networks order) and `Glasser-360`.  See [data](data/) for more details.

## `data`
- The [sub_info](data/sub_info) folder contains the subject information (`sub_id`,`age` and `gender`) used in this study.
- The [fc_variability](data/fc_variability) folder contains the *inter-individual fc variability* matrix for the HCP-D, HCP and ABCC datasets, saved in the `.mat` file.
- The [connectome_matrix](data/connectome_matrix) folder contains the brain networks constructed by *functional connectivity* and *individual variability in structural communicability*.
- The [parcellation_files](data/parcellation_files) folder contains the parcellation files used in this study.

## `functions`
The [functions](functions/) folder contains code and files commonly used in `code`.

## `miscellaneous`

The [miscellaneous](miscellaneous/) folder contains code to batch calculating Spearman's rank correlation between FC variability and other connectomes, along with confidence interval estimation using a bootstrap approach.

## `code`

- The [step_01_inter_individual_fc_variability](step_01_inter_individual_fc_variability/) folder contains codes to construct the individual functional connectivity and estimate the inter-individual variability of functional connectivity. 
- The [step_02_connectional_axis_of_fc_variability](step_02_connectional_axis_of_fc_variability/) folder contains codes to generate results and figures of *Fig. 1. Individual variability in edge-level FC declines along a connectional axis*. 
- The [step_03_white_matter_structural_connectome](step_03_white_matter_structural_connectome/) folder contains codes to generate results and figures of *Fig. 2. Individual variability in structural connectivity communicability is associated with the connectional axis pattern in FC variability across connectome edges*. 
- The [step_04_development_effects](step_04_development_effects/) folder contains codes to generate results and figures of *Fig. 3. Connectional axis of FC variability evolves during youth*.
- The [step_05_cognitive_effects](step_05_cognitive_effects/) folder contains codes to generate results and figures of *Fig. 4. The connectional variability axis pattern is associated with the individual differences in higher-order cognitive functions*.
- The [step_06_psychopathology_effects](step_06_psychopathology_effects/) folder contains codes to generate results and figures of *Fig. 5. The connectional variability axis is linked to both baseline and longitudinal general psychopathology*.

## `wiki`
The detailed description about the codes used in this study, from `step_01_inter_individual_fc_variability` to `step_06_psychopathology_effects`, can be found [here](https://github.com/CuiLabCIBR/Connectional_Variability_axis/wiki).

## `note`
This study primarily used the HCP-D and HCP-YA datasets. The ABCC dataset was only used for the study of psychopathological effects. The current repository does not yet include the code and data used in the supplementary analyses, these will be added in the future. The FC matrices for each individual from the HCP-D and HCP-YA datasets have not been uploaded due to the limitations imposed by the large file sizes.

