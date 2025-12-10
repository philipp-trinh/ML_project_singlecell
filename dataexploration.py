import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Load the single cell data
adata = sc.read_text("in/GSE83523_CE_blastomeres_rawcounts.txt",
                     delimiter="\t", first_column_names=True)

# rename cells
with open("in/GSE83523_series_matrix.txt") as f:
    flines = f.readlines()
lines = [line.split(sep="\t")[1:] for line in flines if
         line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession")]
lines = [(l1.strip('"'), l2.strip('"')) for l1, l2 in zip(*lines)]
lines[-1] = tuple(l.rstrip('\n"') for l in lines[-1])
# otherinfo = [line.split(sep="\t")[1:] for line in flines if line.startswith("!Sample_characteristics_ch1") or line.startswith("")]

# create dictionary with cell names
cell_names = {line[1]: line[0] for line in lines}

# add cell names to adata
adata.var_names = [cell_names[cell] for cell in adata.var_names]

adata.var["embryostage"] = [l.split(sep="_")[0] for l in adata.var_names]
adata.var['embryostage'] = adata.var['embryostage'].astype('category')

adata = adata.T


# load in translation table
trans = pd.read_csv("in/translationfile_celltype.csv", sep="\t")

# plot celltype number
plot = np.unique(sorted(trans["Annotated State"]), return_counts=True)
fig = plt.figure(figsize=(20, 8))
ax = fig.add_subplot(111)
ax.plot(list(map(lambda x: str(x), plot[0])), plot[1])
plt.xlabel("Celltype")
plt.xticks(rotation=90)
plt.ylim(0, 30)
plt.ylabel("Number of cells")
plt.title("Number of cells per celltype")
plt.show()

# boxplot
plt.boxplot(plot[1])
plt.ylabel("Number of cells")
plt.title("Number of cells per celltype")
plt.ylim(0, 30)
plt.show()


sc.pp.neighbors(adata, n_neighbors=20, use_rep='X')
sc.tl.leiden(adata, resolution=0.5)

# PCA
sc.tl.pca(adata, n_comps=50)
sc.pl.pca(adata, color=['embryostage', "leiden"])


sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)
sc.pl.violin(adata.T, ['log1p_mean_counts', 'n_cells_by_counts'], jitter=0.4, multi_panel=True)

# Pre-processing and quality control
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]

# sc.pp.calculate_qc_metrics(adata, inplace=True)
# sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)
# sc.pl.violin(adata.T, ['log1p_mean_counts', 'n_cells_by_counts'], jitter=0.4, multi_panel=True)

plot = np.unique(sorted([int(l.split(sep="_")[0]) for l in adata.obs_names]), return_counts=True)

# lineplot
plt.plot(list(map(lambda x: str(x), plot[0])), plot[1])
plt.xlabel("Embryostage")
plt.xticks(rotation=60)
plt.ylabel("Number of cells")
plt.title("Number of cells per embryostage")
plt.show()

adata = adata[trans["sample name"], :]
obs = list(adata.obs_names)
for i, sample in enumerate(trans["sample name"]):
    for j, s2 in enumerate(adata.obs_names):
        if sample == s2:
            obs[j] = adata.obs_names[j] + "_" + trans["Annotated State"][i].replace(" ", "-")
adata.obs_names = obs
tocsv = pd.DataFrame(adata.X)
tocsv.index = adata.obs_names
tocsv.columns = adata.var_names
tocsv.to_csv("out/transformeddata_onlyannotations.csv")
