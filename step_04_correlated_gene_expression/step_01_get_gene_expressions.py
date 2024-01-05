import abagen
import warnings
from nilearn.datasets import fetch_atlas_schaefer_2018
import pandas as pd

warnings.filterwarnings('ignore', category=FutureWarning)

working_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Gradient/step_04_correlated_gene_expression/'

schaefer = fetch_atlas_schaefer_2018(n_rois=400)
norm_method = 'srs'
expression_donor, report = abagen.get_expression_data(schaefer['maps'], lr_mirror="bidirectional",
                                                      missing="interpolate",
                                                      sample_norm=norm_method, gene_norm=norm_method,
                                                      return_donors=True, return_report=True)

fh = open(working_dir + 'report.txt', 'w', encoding='utf-8')
fh.write(report)
fh.close()

expression_threshold = abagen.correct.keep_stable_genes(list(expression_donor.values()),threshold=0.1, percentile=False)
expression = pd.concat(expression_threshold).groupby('label').mean()
expression.to_csv(working_dir + 'gene_expressions_schaefer400.csv')
