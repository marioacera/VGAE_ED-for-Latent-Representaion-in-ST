import os
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from scanpy.readwrite import read_visium
from scipy.sparse import csr_matrix, hstack

def mk_dir(input_path):
    if not os.path.exists(input_path):
        os.makedirs(input_path)
    return input_path

def adata_preprocess_onthology(i_adata, ont_path, onthology,use_ESM, ESM_column_name):
    if os.path.isfile(ont_path):
        return sc.read(ont_path)
    else: 
        print('Onthology count matrix not found, will be computed')
        if use_ESM:
            gene_name = 'Gene stable ID'
            i_adata.var['Gene name'] = i_adata.var.index
            i_adata.var.index = i_adata.var[ESM_column_name]
            i:adata.var.index.name = None

        else:
            gene_name = 'Gene name'
        data = pd.read_csv(('../data/Onthology/'+ onthology), header = 0)
        print(data.head())
        intersec = len(list(set(i_adata.var.index).intersection(list(data.drop_duplicates(subset = gene_name)[gene_name]))))
        if intersec == 0:
            print('Onthology reference not matching sample genes')
            exit()
        print('Number of genes in the sample found in the onthology reference:',intersec)

        ont = sc.AnnData(obs = i_adata.obs, var = pd.DataFrame(data.drop_duplicates(subset = 'GO term name')['GO term name']))
        ont.var.index = ont.var['GO term name']
        ont.uns['spatial'] = i_adata.uns['spatial']
        ont.obsm['spatial'] = i_adata.obsm['spatial']
        ont.var['positions'] = np.arange(0,len(ont.var.index.to_list()))
        i_adata.var['positions'] = np.arange(0,len(i_adata.var.index.to_list()))
        ont.X = np.zeros((len(i_adata.obs.index.to_list()),len(ont.var.index.to_list())))
        print('                                                                         ')
        for gene in i_adata.var.index:
            i = i_adata.var.loc[gene].positions
            if i%1000 == 0: print('\n',"\033[A                             \033[A",'Progress :',i/len(i_adata.var.index.to_list()))
            aux_data = data[data[gene_name] == gene]
            try:
                counts_to_sum = i_adata.X.todense()[:,i:i+1]
            except:
                counts_to_sum = i_adata.X[:,i:i+1]

            aux_data = aux_data[pd.isna(aux_data['GO term name'])== False]
            for GO_term in aux_data['GO term name']:
                j = ont.var.loc[GO_term].positions 
                ont.X[:,j:j+1] += counts_to_sum
        ont.var['GO termn'] = ont.var.index
        ont.var.index = ont.var.positions.map(str)
        ont.var.index.name = None
        ont.write(ont_path, compression = 'gzip')
        return ont


def adata_preprocess(i_adata, min_cells=3, pca_n_comps=300, use_onthology = False, ont_count_data_name = None , onthology_type = 'Human_Onthology.txt', use_ESM = 'False', ESM_column_name = 'gene_ids'):
    print('===== Preprocessing Data ')
    sc.pp.filter_genes(i_adata, min_cells=min_cells)
    adata_X = sc.pp.normalize_total(i_adata, target_sum=1, exclude_highly_expressed=True, inplace=False)['X']
    adata_X = sc.pp.scale(adata_X)
    if use_onthology == True:

        ont = adata_preprocess_onthology(i_adata, ont_count_data_name, onthology_type, use_ESM, ESM_column_name)
        sc.pp.filter_genes(ont, min_cells = min_cells)
        ont.var_names_make_unique()
        ont_X = sc.pp.normalize_total(ont, target_sum=1, exclude_highly_expressed=True, inplace=False)['X']
        ont_X = sc.pp.scale(ont_X,)
        ont_X = sc.pp.pca(ont_X, n_comps=int(pca_n_comps/2)) 
        adata_X = sc.pp.pca(adata_X, n_comps=int(pca_n_comps/2))
        return np.concatenate((np.array(adata_X),np.array(ont_X)), axis = 1)
    else:
        adata_X = sc.pp.pca(adata_X, n_comps=pca_n_comps) 
    #sc.pp.pca(i_adata, n_comps=300)
    #sc.pp.highly_variable_genes(i_adata, n_top_genes= pca_n_comps)
        return adata_X
    #i_adata[:,i_adata.var.highly_variable]
    
def adata_preprocess_ont(i_adata, min_cells=3, pca_n_comps=300, use_onthology = False, ont_count_data_name = None , onthology_type = 'Human_Onthology.txt', use_ESM = 'False', ESM_column_name = 'gene_ids'):
    print('===== Preprocessing Data ')
    ont = adata_preprocess_onthology(i_adata, ont_count_data_name, onthology_type, use_ESM, ESM_column_name)
    sc.pp.filter_genes(ont, min_cells = min_cells)
    ont.var_names_make_unique()
    ont_X = sc.pp.normalize_total(ont, target_sum=1, exclude_highly_expressed=True, inplace=False)['X']
    ont_X = sc.pp.scale(ont_X,)
    ont_X = sc.pp.pca(ont_X, n_comps=int(pca_n_comps)) 

    #sc.pp.pca(i_adata, n_comps=300)
    #sc.pp.highly_variable_genes(i_adata, n_top_genes= pca_n_comps)
    return ont_X
    #i_adata[:,i_adata.var.highly_variable]
   
def adata_preprocess_ont_hvg(i_adata, min_cells=3, pca_n_comps=300, use_onthology = False, ont_count_data_name = None , hvg = 2000, onthology_type = 'Human_Onthology.txt', use_ESM = 'False', ESM_column_name = 'gene_ids'):
    print('===== Preprocessing Data ')
    ont = adata_preprocess_onthology(i_adata, ont_count_data_name, onthology_type, use_ESM, ESM_column_name)
    ont = adata_preprocess_onthology(i_adata, ont_count_data_name, onthology_type, use_ESM, ESM_column_name)
    sc.pp.filter_genes(ont, min_cells=min_cells)
    sc.pp.normalize_total(ont, target_sum=1, exclude_highly_expressed=True)
    sc.pp.log1p(ont)
    sc.pp.highly_variable_genes(ont,  n_top_genes = hvg)
    ont = ont[:, ont.var.highly_variable]
    sc.pp.scale(ont)
    ont_X = ont.X
    ont_X = sc.pp.pca(ont_X, n_comps=int(pca_n_comps)) 
    #sc.pp.pca(i_adata, n_comps=300)
    #sc.pp.highly_variable_genes(i_adata, n_top_genes= pca_n_comps)
    return ont_X

def adata_preprocess_hvg(i_adata, min_cells=3, pca_n_comps=300, use_onthology = False, ont_count_data_name = None, hvg = 2000, onthology_type = 'Human_onthology.txt',use_ESM = 'False', ESM_column_name = 'gene_ids'):
    print('===== Preprocessing Data ')
    sc.pp.filter_genes(i_adata, min_cells=min_cells)
    sc.pp.normalize_total(i_adata, target_sum=1, exclude_highly_expressed=True)
    sc.pp.log1p(i_adata)
    sc.pp.highly_variable_genes(i_adata, n_top_genes = hvg)
    i_adata = i_adata[:, i_adata.var.highly_variable]
    sc.pp.scale(i_adata)
    if use_onthology == True:
        ont = adata_preprocess_onthology(i_adata, ont_count_data_name, onthology_type, use_ESM, ESM_column_name)
        sc.pp.filter_genes(ont, min_cells=min_cells)
        sc.pp.normalize_total(ont, target_sum=1, exclude_highly_expressed=True)
        sc.pp.log1p(ont)
        sc.pp.highly_variable_genes(ont,  n_top_genes = hvg)
        ont = ont[:, ont.var.highly_variable]
        sc.pp.scale(ont)
        ont_X = ont.X
        adata_X = i_adata.X
        ont_X = sc.pp.pca(ont_X, n_comps=int(pca_n_comps/2)) 
        adata_X = sc.pp.pca(adata_X, n_comps=int(pca_n_comps/2))
        return np.concatenate((np.array(adata_X),np.array(ont_X)), axis = 1)
    else:
        adata_X = i_adata.X.copy()
        adata_X = sc.pp.pca(adata_X, n_comps=pca_n_comps)

    #sc.pp.pca(i_adata, n_comps=300)
    #sc.pp.highly_variable_genes(i_adata, n_top_genes= pca_n_comps)
        return adata_X




def load_ST_file(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True, file_Adj=None):
    adata_h5 = sc.read_visium(file_fold, load_images=load_images, count_file=count_file)
    adata_h5.var_names_make_unique()

    if load_images is False:
        if file_Adj is None:
            file_Adj = os.path.join(file_fold, "spatial/tissue_positions_list.csv")

        positions = pd.read_csv(file_Adj, header=None)
        positions.columns = [
            'barcode',
            'in_tissue',
            'array_row',
            'array_col',
            'pxl_col_in_fullres',
            'pxl_row_in_fullres',
        ]
        positions.index = positions['barcode']
        adata_h5.obs = adata_h5.obs.join(positions, how="left")
        adata_h5.obsm['spatial'] = adata_h5.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()
        adata_h5.obs.drop(columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True)

    print('adata: (' + str(adata_h5.shape[0]) + ', ' + str(adata_h5.shape[1]) + ')')
    return adata_h5
