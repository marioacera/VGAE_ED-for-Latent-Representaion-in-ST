import os
import torch
import argparse
import warnings
import numpy as np
import pandas as pd
from graph_utils import graph_construction
from utils import mk_dir, adata_preprocess_hvg, load_ST_file,adata_preprocess, adata_preprocess_ont,adata_preprocess_ont_hvg
import anndata
from train import Train
from sklearn import metrics
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

warnings.filterwarnings('ignore')
torch.cuda.cudnn_enabled = False
np.random.seed(0)
torch.manual_seed(0)
torch.cuda.manual_seed(0)
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
print('===== Using device: ' + device)

# ################ Parameter setting
parser = argparse.ArgumentParser()

#
parser.add_argument('--use_onthology', type=bool, default=False, help='Using onthology, reference based, profiling')
parser.add_argument('--use_ESM', type=bool, default=True, help='When computing pathways count matrix, use ESM name instead of gene symbol')
parser.add_argument('--ESM_column_name', type=str, default='gene_ids', help='Column where the ESM name is located')
parser.add_argument('--use_only_hvg', type=bool, default=False, help='Using highly variable features')
parser.add_argument('--n_hvg', type=int, default=5000, help='Nº of highly variable features')

#___________________ Model Parameters ______________________
parser.add_argument('--k', type=int, default=15, help='parameter k in spatial graph')
parser.add_argument('--knn_distanceType', type=str, default='euclidean', help='graph distance type: euclidean/cosine/correlation')
parser.add_argument('--epochs', type=int, default=300, help='Number of epochs to train.')
parser.add_argument('--cell_feat_dim', type=int, default=300, help='Dim of PCA')
parser.add_argument('--feat_hidden1', type=int, default=100, help='Dim of DNN hidden 1-layer.')
parser.add_argument('--feat_hidden2', type=int, default=20, help='Dim of DNN hidden 2-layer.')
parser.add_argument('--gcn_hidden1', type=int, default=32, help='Dim of GCN hidden 1-layer.')
parser.add_argument('--gcn_hidden2', type=int, default=8, help='Dim of GCN hidden 2-layer.')
parser.add_argument('--p_drop', type=float, default=0.2, help='Dropout rate.')
parser.add_argument('--using_dec', type=bool, default=True, help='Using DEC loss.')
parser.add_argument('--using_mask', type=bool, default=False, help='Using mask for multi-dataset.')
parser.add_argument('--feat_w', type=float, default=10, help='Weight of DNN loss.')
parser.add_argument('--gcn_w', type=float, default=0.1, help='Weight of GCN loss.')
parser.add_argument('--dec_kl_w', type=float, default=1, help='Weight of DEC loss.')
parser.add_argument('--gcn_lr', type=float, default=0.01, help='Initial GNN learning rate.')
parser.add_argument('--gcn_decay', type=float, default=0.01, help='Initial decay rate.')

parser.add_argument('--dec_cluster_n', type=int, default=8, help='DEC cluster number.')
parser.add_argument('--dec_interval', type=int, default=20, help='DEC interval nnumber.')
parser.add_argument('--dec_tol', type=float, default=0.00, help='DEC tol.')

# ______________ Eval clustering Setting _________
parser.add_argument('--eval_resolution', type=int, default=0.5, help='Eval cluster number.')
parser.add_argument('--eval_graph_n', type=int, default=15, help='Eval graph kN tol.') 


params = parser.parse_args()
params.device = device

# ################ Path setting
data_root = '/home/marioam/Documents/Deep Learning Applications/Autoencoder-Onthology/data/'
data_name = 'Human_brain_151673'
ont_count_data_name = '151673_ont.h5ad'
onthology_type = 'Human_Onthology.txt' # Files allocated in ../data/Onthology/ #### FILES DOWNLOADED FROM: http://www.ensembl.org/biomart/martview/83827fd6fb9a38932799e5d2b96c4685 
save_fold = os.path.join('../results/', data_name)

# ################## Load data
adata_h5 = load_ST_file(file_fold=os.path.join(data_root, data_name))
adata_h5.var_names_make_unique()
if params.use_only_hvg:
    adata_X_ = adata_preprocess_ont_hvg(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim, use_onthology = params.use_onthology, ont_count_data_name = ont_count_data_name, hvg = params.n_hvg, onthology_type = onthology_type , use_ESM = params.use_ESM, ESM_column_name = params.ESM_column_name)
else:
    adata_X_ = adata_preprocess_ont(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim, use_onthology = params.use_onthology, ont_count_data_name = ont_count_data_name, onthology_type = onthology_type, use_ESM = params.use_ESM, ESM_column_name = params.ESM_column_name)

graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)
params.cell_num = adata_h5.shape[0]
params.save_path = mk_dir(save_fold)
print('==== Graph Construction Finished')
print(adata_X_.shape)
#adata_X_ = adata_X_[:,0:params.cell_feat_dim]
#adata_X = adata_X_.X
# ################## Model training
sed_net = Train(adata_X_, graph_dict, params)
if params.using_dec:
    sed_net.train_with_dec()
else:
    sed_net.train_without_dec()
sed_feat, de_feat,_, feat_x, _ = sed_net.process()

np.savez(os.path.join(params.save_path, "VGAE_result.npz"), sed_feat=sed_feat, deep_Dim=params.feat_hidden2)

# ################## Result plot
adata_sed = adata_h5
#adata_sed.X = de_feat
adata_sed.uns['spatial'] = adata_h5.uns['spatial']
adata_sed.obsm['spatial'] = adata_h5.obsm['spatial']

adata_sed.obsm['red_X'] = sed_feat
#adata_X_.obsm['red_X'] = sed_feat

sc.pp.neighbors(adata_sed, n_neighbors=params.eval_graph_n, use_rep = 'red_X')
sc.tl.umap(adata_sed)
sc.tl.leiden(adata_sed, key_added="VGAE_leiden", resolution=params.eval_resolution)
sc.pl.umap(adata_sed, color = ['VGAE_leiden'],  save = '_Human_brain_leiden_regions_ONT_05resolution.png')
sc.pl.spatial(adata_sed, img_key="hires", color=['VGAE_leiden'], save = '_Human_brain_leiden_regions_ONT_05esolution.png')
#sc.pl.spatial(adata_sed, img_key="hires", color=[adata_sed.var.index[0],adata_sed.var.index[1]], save = "Gene_patterns_refined_plot.png")
#plt.savefig(os.path.join(params.save_path, "Gene_patterns_refined_plot.png"), bbox_inches='tight', dpi=150)

#adata_X_.obsm['X_umap'] = adata_sed.obsm['X_umap']
#sc.pl.spatial(adata_X_, img_key="hires", color=[adata_X_.var.index[0],adata_X_.var.index[1]], save = "Gene_patterns_reference_plot.png")
#plt.savefig(os.path.join(params.save_path, "Gene_patterns_reference_plot.png"), bbox_inches='tight', dpi=150)


df_result = pd.DataFrame(adata_sed.obs['VGAE_leiden'], columns=['VGAE_leiden'])
df_result.to_csv(os.path.join(params.save_path, "VGAE_leiden_n_"+str(params.eval_resolution)+"_result.tsv"),
                 sep='\t', index=False)
