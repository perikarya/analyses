pip install 'scanpy[leiden]' --quiet

pip install 'gseapy' --quiet

pip install 'seaborn' --quiet

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from matplotlib import cm
import numpy as np
import gseapy
import seaborn as sns

dat = sc.read_h5ad("/content/drive/MyDrive/colab_files/muris/tabula-muris-senis-facs-processed-official-annotations-Brain_Non-Myeloid.h5ad")

print(dat)
print(dat.shape)
print(dat.raw.shape)

logged = sc.pp.log1p(dat, base=2, copy=True)

logged.raw = logged
print(logged)
print(logged.shape)
print(logged.raw.shape)

sc.pl.umap(logged, color='free_annotation')
sc.pl.tsne(logged, color='free_annotation')

sc.pl.umap(logged, color='cell_ontology_id')
sc.pl.tsne(logged, color='cell_ontology_id')

sc.pl.umap(logged, color='cell_ontology_class', palette=sns.color_palette("husl", 17))
sc.pl.tsne(logged, color='cell_ontology_class')

sc.tl.rank_genes_groups(logged, 'louvain', method='t-test_overestim_var', key_added = "t-test_louvain")
sc.pl.rank_genes_groups(logged, n_genes=25, sharey=False, key = "t-test_louvain")

sc.tl.rank_genes_groups(logged, 'leiden', method='t-test', key_added = "t-test_leiden")
sc.pl.rank_genes_groups(logged, n_genes=25, sharey=False, key = "t-test_leiden")

sc.tl.rank_genes_groups(logged, 'cell_ontology_class', method='t-test', key_added = "t-test_cell_ontology_class")
sc.pl.rank_genes_groups(logged, n_genes=25, sharey=False, key = "t-test_cell_ontology_class")
