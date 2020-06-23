import DukeDS
import os

quarter = os.getenv('DDH_QUARTER', '20Q1')

ddh_files = [
  'gene_summary.Rds',
  '{}_achilles.Rds'.format(quarter),
  '{}_expression_join.Rds'.format(quarter),
  '{}_subcell.Rds'.format(quarter),
  'pathways.Rds',
  'surprise_genes.Rds',
  '{}_pubmed_concept_pairs.Rds'.format(quarter),
  '{}_na_cutoff.Rds'.format(quarter),
  'sd_threshold.Rds',
  'achilles_lower.Rds',
  'achilles_upper.Rds',
  'mean_virtual_achilles.Rds',
  'sd_virtual_achilles.Rds',
  'master_bottom_table.Rds',
  'master_top_table.Rds',
  'master_positive.Rds',
  'master_negative.Rds',
  ]

for ddh_file in ddh_files:
  print(ddh_file)
  DukeDS.download_file('ddh-data', 'data/{}'.format(ddh_file))
