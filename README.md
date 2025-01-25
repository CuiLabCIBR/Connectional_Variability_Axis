# Connectional_Variability_Axis
Data and codes for our paper **"Connectional Axis of Individual Functional Variability: Patterns, structural correlates, and implications for development and cognition**.

*Inter-individual fc variability* matrices for the Human Connectome Project (HCP)-development  (**HCP-D**),  unrelated HCP-young adult (**HCP-YA**, typically referred as **'HCP'**) and Youth Executive function and Neurodevelopment (**YEN**) datasets, and the connectome matrices of *functional connectivity* and *individual variability in structural communicability* are obtained using `Schaefer-400` (7 Networks order) and `Glasser-360`.  See [data](data/) for more details.

## `data`
- The [sub_info](data/sub_info) folder contains the subject information (`sub_id`,`age` and `sex`) used in this study. 
- The [fc_variability](data/fc_variability) folder contains the *inter-individual fc variability* matrix for the HCP-D, HCP and YEN datasets, saved in the `.mat` file.
- The [connectome_matrix](data/connectome_matrix) folder contains the brain networks constructed by *functional connectivity* and *individual variability in structural communicability*.
- The [parcellation_files](data/parcellation_files) folder contains the parcellation files used in this study.

## `functions`
The [functions](functions/) folder contains code and files commonly used in `code`.

## `code`

- The [step_01_inter_individual_fc_variability](step_01_inter_individual_fc_variability/) folder contains codes to construct the individual functional connectivity and estimate the inter-individual variability of functional connectivity. 
- The [step_02_connectional_axis_of_fc_variability](step_02_connectional_axis_of_fc_variability/) folder contains codes to generate results and figures of *Fig. 1. Individual variability in edge-level FC declines along a connectional axis*. 
- The [step_03_structural_connectome_variability](step_03_structural_connectome_variability/) folder contains codes to generate results and figures of *Fig. 2. Individual variability in structural connectivity communicability is associated with the connectional axis pattern in FC variability across connectome edges*. 
- The [step_04_developmental_effects](step_04_developmental_effects/) folder contains codes to generate results and figures of *Fig. 3. Connectional axis of FC variability evolves during youth*.
- The [step_05_cognitive_effects](step_05_cognitive_effects/) folder contains codes to generate results and figures of *Fig. 4. The connectional variability axis pattern is associated with the individual differences in higher-order cognitive functions*.
- The scripts to generate the *Fig. 5. Replication in an independent developmental dataset*. can be found in the sensitivity analyses within the aforementioned folders from `step_01_inter_individual_fc_variability` to `step_05_cognitive_effects`.

## `wiki`
The detailed description about the codes used in this study, from `step_01_inter_individual_fc_variability` to `step_05_cognitive_effects`, can be found [here](https://github.com/CuiLabCIBR/Connectional_Variability_axis/wiki).

