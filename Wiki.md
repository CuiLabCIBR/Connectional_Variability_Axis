<a name="top"></a> Welcome to the Connectional_Variability_Axis wiki!

## Abstract
The human cerebral cortex exhibits intricate interareal functional synchronization at the macroscale, with substantial individual variability in these functional connections. However, the spatial organization of functional connectivity (FC) variability across the human connectome edges and its significance in cognitive development remain unclear. Here, we identified a connectional axis in the edge-level FC variability. The variability declined continuously along this axis from within-network to between-network connections, and from the edges linking association networks to those linking the sensorimotor and association networks. This connectional axis of functional variability is associated with spatial pattern of structural connectivity variability. Moreover, the connectional variability axis evolves in youth with an increasing flatter axis slope. We also observed that the slope of connectional variability axis was positively related to the performance in the higher-order cognition. Together, our results reveal a connectional axis in functional variability that is linked with structural connectome variability, refines during development, and is relevant to cognition.

## Run the codes in serieal
### step_01_inter_individual_fc_variability
This folder contains codes to extract the parcel-level BOLD signal based on the `Schaefer-400` atlas, calculate the FC, and estimate the individual variability of FC. These codes are almost identical for HCP-D and HCP (i.e., HCP-YA), except for differences in some parameters, such as the session list. We estimated the inter-individual variability and intra-individual variability based on a linear fixed-effects (LME) model, and more details are described in [FC variability](#fc-variability). 

> **Notes:** (As the original FC matrices from all participants are too large and cannot be uploaded to GitHub, so if you want to replicate this project, you can skip step_01_individual_fc_variability and directly use our FC variability results, and start from step_02_connectional_axis_of_fc_variability.)

(1) step_01_calculate_fc_hcp/hcpd.m
> This script is used to calculate the FC based on the `Schaefer-400` atlas. The HCP minimal preprocessed fMRI data were post-processed using **xcp-d** (https://xcp-d.readthedocs.io/). Post-processed steps include nuisance regression, filtering and smooth. Post-processed data `*_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii` for each participant was used, and then we sorted the 400 cortical regions based on Yeo_7Networks. 
>
> For HCP-D, **Four rest runs** were concatenated and separated into **eight runs** to estimate `inter-individual variability` and `intra-individual variability`, respectively. For HCP, four rest runs were concatenated and separated into **twelve runs** to estimate `inter-individual variability` and `intra-individual variability`, respectively.
>
> Finally, the FC was calculated as the Pearson correlation between timeseries from region pairs. `Schaefer-400` atlas was used, generating a `400×400 FC matrix` for each participant, and there are **79,800** unique edges (399×400/2) in each matrix. The r values were converted to z values using Fisher r-to-z transformation.

(2) step_02_load_fc_hcp/hcpd.m

> This script is used to load the FC from all participants and save the results in a 3D matrix, dim: participant × session × edge.

(3) step_03_prepare_lme_data.m

> This script is used to prepare data for fc variability calculation based on the LME model in R. The data is saved as a 2D matrix, dim: (session × participant) × edge. Taking the HCP-D dataset as an example, there are 486 participants, each participant has 8 runs, and each run generates a FC matrix with 79,800 unique edges. Therefore, the dimension of the data matrix is 3,888×79,800. This script also generates the `subID` and `seesion` used in the LME model.

(4) step_04_lme.R

> This script is used to estimate the `inter-individual variability` and `intra-individual variability` of FC edges using a linear mixed-effects (LME) model.

(5) step_05_get_inter_individual_fc_variability.m

> The script is used to get the inter-individual FC variability matrix based on the `Schaefer-400` atlas.

---
### step_02_connectional_axis_of_fc_variability
This folder contains codes to generate figures and results of `Figure 1. Individual variability in edge-level FC declines along a connectional axis`. 

Parcels from the`Schaefer-400` were mapped onto seven canonical large-scale functional networks from the Yeo atlas. We excluded the limbic network in the following analyses, as previous studies have consistently reported substantial signal loss in this network, especially in the orbitofrontal and ventral temporal cortices. Our analysis contained 374 parcels from six functional networks, including the visual (VS), somatomotor (SM), dorsal attention (DA), ventral attention (VA), frontoparietal (FP), and default mode (DM) networks.

The connectional axis of this variability matrix can be uniquely represented by six within-network FC variability values (diagonal elements) and 15 between-network FC variability values (nondiagonal elements). Twenty-one unique FC variability values were defined as the profile of connectional variability axis.

(1) step_01_connectional_axis_of_fc_variability.m
> We plotted the inter-individual FC variability matrix at the nodal-level and network-level for both the HCP-D and HCP datasets. The similarity between the two datasets was measured by Spearman's rank correlation. We grouped all edges into three types, including within-sensorimotor (S-S), within-association (A-A), and between association-sensorimotor (A-S), and compared their mean differences in FC variability using a permutation test (10,000 iterations).

(2) step_02_ci_bootstrap.R

> The confidence interval of Spearman's rank correlation was estimated using a bootstrap approach. 

(3) step_05_barplot_connectional_variability_axis_single_net.R
> Bar plot for the connectional variability axis represented by within-network to between-network FC variability for each functional network (e.g., DM-DM, DM-FP, DM-VA, DM-DA, DM-SM, and DM-VS).

(4) step_06_barplot_connectional_variability_axis_all_net.R

> Bar plot for the connectional variability axis represented by all 21 within-network to between-network FC variability for each functional network (e.g., DM-DM, DM-FP, DM-VA, DM-DA, DM-SM, and DM-VS).

(5) step_07_boxplot_S_A_axis.R

> Box plot for the FC variability values within-sensorimotor (S-S), within-association (A-A), and between sensorimotor-association (S-A).

---
### step_03_structural_connectome_variability
This folder contains codes to generate figures and results of `Figure 2. Individual variability in structural connectivity communicability is associated with the connectional axis pattern in FC variability across connectome edges`. The probabilistic tractography was performed using MRtrix3, with 10 million streamlines, length range from 30mm to 250mm, FOD power = 0.33. Next, a structural connectivity matrix was constructed using *tck2connectome* with the `Schaefer-400` atlas. The weight of each structural connectivity indicates the number of streamlines connecting two regions. Edge weights were normalized by dividing the average volume of the two regions, and then log-transformed. Particularly, edges with a weight of zero for any participant were set to zero for all participants.

(1) step_01_structural_communicability_variability.m
> Considering that the white matter structural connectome accounts only for **direct** functional communication, we calculated the communicability of the structural connectivity matrix for each participant to include both the direct and **indirect** pathways. `Communicability` quantifies the capacity of two regions to communicate through other regions by pathways of all possible topological lengths. Specifically, let *A* denotes the weighted adjacency matrix of the individual structural connectome. The communicability between two regions *i* and *j* is calculated as $$(\exp(D^{-1/2}A D^{-1/2}))_{i j}$$, where $$D\,=\,d i a g(\Sigma_{k=1}^{N}\,a_{i k})$$ is a diagonal matrix in which the diagonal elements represent the strength of each node. We calculated the `mean absolute deviation` of structural communicability across all participants for each edge and defined it as the **inter-individual variability of structural communicability** across all participants.
>
> Then, we evaluated the alignment between FC variability and structural communicability variability across all edges by computing Spearman’s rank correlation of the connection strength between the communicability variability and FC variability.

(2) step_02_ci_bootstrap.R

> The confidence interval of Spearman's rank correlation was estimated using a bootstrap approach. 

(3) step_03_density_plot.py
> Density plot for the alignment between `FC variability` and `Structural communicability variablity` for both HCP-D and HCP at the nodal level. 

---
### step_04_developmental_effects
This folder contains codes to generate figures and results of `Figure 3. Connectional axis of FC variability evolves during youth`.

(1) step_01_development_connectional_variability_axis.m
> To explore the maturation of connectional variability axis in youth, we divided the HCP-D participants (aged 8–21 years old) into multiple subgroups and calculated the inter-individual variability matrix for each group. We balanced the number of groups and group sizes using a **sliding window** approach to ensure sufficient statistical power and accurate estimation of individual variability. Specifically, 486 participants from the HCP-D dataset were sorted in ascending order by age and then divided into overlapping groups. We set the window length to 50 participants with a step size of five participants, and obtained 88 groups.
>
> We evaluated how the connectional axis of individual FC variability developed in youth. The connectional variability axis represents the ranking of the 21 within-network and between-network average variabilities. We calculated the slope of the connectional variability axis using linear regression to quantify the heterogeneity of the axis, and estimated the `axis slope` for each of the 88 groups. Specifically, the linear model was constructed as Y = β0 + β1*X, where Y represents the 21 ranked FC variability values, and the X is the rank ranging from 1 to 21 (scaled from 0 to 1 by dividing each rank by 21), β0 is the intercept, and β1 measures the slope of the connectional variability axis.

(2) step_02_plot_gam_development.R
> To model both linear and nonlinear developmental effects of  `axis slope`, we used a generalized additive model (GAM) with penalized splines, which estimates nonlinearities using restricted maximum likelihood (REML) and penalizes nonlinearity to avoid overfitting data. GAM was fit with the variability axis alignment was modeled as the dependent variable. The group age was modeled as a smooth term while group sex and group in-scanner head motion as model covariates as follows:
>
> $$\text {Variability axis slope} \sim s(\text {Age, } k=3)+\text {Sex}+\text {Motion}$$
>
> where $s()$ is the spline basis function, and $k$ determines the maximum basis complexity. As the analyses were performed at the group level, a group average age was assigned for each group. The sex of each group was computed as the ratio of males to females. The head motion was calculated by averaging the mean FD across all individuals and all runs within each group.
>
> The significance of the association between `axis slope` and age was assessed through analysis of variance (ANOVA) that compared the full GAM model to a nested, reduced model with no age term. A significant result indicated that the residual deviance was significantly lower when a smooth term for age was included in the model, as assessed using the chi-squared test statistic. To determine the overall magnitude and direction of the association between axis alignment and age, we calculated the partial *R*<sup>2</sup> between the full GAM and reduced models (effect magnitude), and signed the partial *R*<sup>2</sup> by the sign of the age coefficient from an equivalent linear model (effect direction).

(3) step_03_gam_development_single_net_edge.R

> To determine which connections drive the change in axis slope, we constructed the same GAM model for each network-level FC and estimated the developmental effects of its FC variability.

(4) step_04_plot_edge_fc_var_age_effects.m

> Matrix plot for the developmental effects of the 21 network-level FC variability values.

(5) step_05_barplot_age_effects_axis_all_net.R

> Bar plot for the developmental effects of the 21 network-level FC variability values.

---
### step_05_cognitive_effects

This folder contains codes to generate figures and results of `Figure 4. The connectional variability axis pattern is associated with the individual differences in higher-order cognitive functions`. We used the `composite score of fluid cognition` from the NIH Toolbox Cognition Battery to quantify the participants’ cognitive ability. The fluid cognition composite score was obtained by averaging the normalized scores from multiple cognitive tasks, including flanker inhibition, dimensional change card sort (flexibility), picture sequence memory, list sorting working memory, and pattern comparison. This analysis included 352 participants from the HCP-D and 274 from the HCP-YA with fluid cognition composite scores. Taking the HCP-D as examples:

(1) step_01_prepare_data_cognition_hcp/hcpd.m

> We ranked the participants based on their cognitive scores and then employed a sliding window approach (window length = 50, step size = 5) to divide them into multiple groups with increasing cognitive scores. In total, 61 groups were obtained for the HCP-D dataset.

(2) step_02_plot_gam_cognition.R

> A GAM analysis was performed to evaluate the relationship between variability axis slope and fluid cognition composite scores while controlling the linear and nonlinear effects of age, as well as the effects of sex and in-scanner motion. The equation for GAM is as follows:
>
> $$\text {Variability axis slope} \sim \text {Cognition}+s(\text {Age, } k=3)+\text {Sex}+\text {Motion}$$
>
> Similar to the developmental analyses, the significance of the association between the axis slope and cognitive performance was assessed using an ANOVA that compared the full GAM model to a nested, reduced model with no cognition term. To evaluate the effect size, we calculated the partial *R*<sup>2</sup> between the full GAM model and the reduced model and signed the partial *R*<sup>2</sup> using the sign of the cognition coefficient (T-value) from the full GAM model.

(3) step_03_gam_cognition_single_net_edge_hcp/hcpd.R

> To determine which connections drive the change in axis slope, we constructed the same GAM model for each network-level FC and estimated the cognitive effects of its FC variability.

(4) step_04_plot_edge_fc_var_cog_effects.m

> Matrix plot for the cognitive effects of the 21 network-level FC variability values.

(5) step_05_barplot_cog_effects_axis_all_net.R

> Bar plot for the cognitive effects of the 21 network-level FC variability values.

---

## FC variability

In this study, `inter-individual variability` and `intra-individual variability` of each FC edge (69,751 unique edges in total) were calculated using a **linear mixed-effects (LME)** model for the HCP-YA and HCP-D datasets36. The LME model was implemented via R package Rex (https://github.com/TingsterX/Reliability_Explorer), which was designed to facilitate the robust measurements of individual variations. This model captures fixed and random effects, assuming that the observed response variable and residual term follow a normal distribution with zero mean and a specific variance. Specifically, for each FC edge, the LME model is expressed as:
$$FC_{i t}=\mu_{0}+\lambda_{i}+\alpha_{t}+\epsilon_{i t}, \text { where } \lambda_{i} \sim N\left(0, \sigma_{\lambda}^{2}\right), \epsilon_{t} \sim N\left(0, \sigma_{\epsilon}^{2}\right) \text {.}$$
Here, $i$ identifies the individual, *t* indicates the session, $\lambda_{i}$ and $\alpha_{t}$ measures random effect and fixed effect, respectively, and $\epsilon_{i t}$ is the residual term. The observed individual variation $\sigma_{FC}^{2}$ between $FC_{i t}$ can be decomposed into real inter-individual variation $\sigma_{b}^{2}$ across participants and the intra-individual variations $\sigma_{w}^{2}$ across sessions (captured by the residual variation $\sigma_{\epsilon}^{2}$). In this work, we used intra-class correlation (ICC, $\frac{\sigma_{b}^{2}}{\sigma_{w}^{2}+\sigma_{b}^{2}}$) to measure the adjusted inter-individual variability while accounting for intra-individual variability.

## Null models

To evaluate the significance of the correlation between the two connectomes across edges, we performed permutation testing by generating null connectomes, as in the spatial permutation test, to compare cortical properties. The spin test was designed to assess the spatial similarity between two cortical maps while considering autocorrelations in the brain data. This is achieved by projecting the brain surface onto a spherical space, rotating the sphere, and projecting it back onto the surface. Here, we utilized the Schaefer atlas, comprising 400 cortical regions with corresponding node IDs ranging from 1 to 400 (i.e., [1, 2, ..., 400]). To create the null connectomes, we shuffled the 400 cortical regions by rotating the Schaefer atlas. We then reordered the rows and columns of the original matrix based on the shuffled node IDs (e.g., [134, 25, ..., 365]). This procedure preserved the mean and variance of the matrix while ensuring that the regional profile (i.e., the matrix column, 400×1) included all 400 cortical regions. 

Repeating this process, we constructed a null distribution by generating 10,000 randomized individual variability connectomes. Next, we calculated Spearman’s rank correlations between the randomized variability and other neurobiological connectomes, including the connectomes of structural communicability variability, correlated gene expression, and correlated neurotransmitter receptor expression. We compared the correlation values obtained from the empirical individual variability matrix with those acquired using the null connectomes to determine the significance level, and the *P* value of this permutation testing (*P<sub>Perm</sub>*) was reported.

## Notes

> The scripts for sensitivity analyses and replication using the YEN dataset can be found within each folder described above.

[Back to the top.](#top)