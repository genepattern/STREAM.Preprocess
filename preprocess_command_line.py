#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
warnings.filterwarnings('ignore')

__tool_name__='STREAM'
print('''
   _____ _______ _____  ______          __  __ 
  / ____|__   __|  __ \|  ____|   /\   |  \/  |
 | (___    | |  | |__) | |__     /  \  | \  / |
  \___ \   | |  |  _  /|  __|   / /\ \ | |\/| |
  ____) |  | |  | | \ \| |____ / ____ \| |  | |
 |_____/   |_|  |_|  \_\______/_/    \_\_|  |_|
                                               
''',flush=True)

import stream as st
import argparse
import multiprocessing
import os
from slugify import slugify
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import sys
mpl.use('Agg')
mpl.rc('pdf', fonttype=42)

os.environ['KMP_DUPLICATE_LIB_OK']='True'


print('- STREAM Single-cell Trajectory Reconstruction And Mapping -',flush=True)
print('Version %s\n' % st.__version__,flush=True)
    

def main():
    sns.set_style('white')
    sns.set_context('poster')
    parser = argparse.ArgumentParser(description='%s Parameters' % __tool_name__ ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", "--matrix", dest="input_filename",default = None, help="input file name", metavar="FILE")
    parser.add_argument("-l", "--cell_labels",dest="cell_label_filename", default=None,help="filename of cell labels")
    parser.add_argument("-c","--cell_labels_colors",dest="cell_label_color_filename", default=None, help="filename of cell label colors")
    parser.add_argument("--log2",dest="flag_log2", action="store_true", help="perform log2 transformation")
    parser.add_argument("--norm",dest="flag_norm", action="store_true", help="normalize data based on library size")
    parser.add_argument("-o","--output_folder",dest="output_folder", default=None,   help="Output folder")
    parser.add_argument("-rmt","--remove_mt_genes",dest="flag_remove_mt_genes", action="store_true", default=False, help="Remove Mitochondrial genes")
    parser.add_argument("-mcg","--min_count_genes",dest="min_count_genes", type = int, default=None,  help="filter cells with less than this many genes")
    parser.add_argument("-mpg","--min_percent_genes",dest="min_percent_genes", type=float, default=None,  help="The minimum percent genes")
    parser.add_argument("-mpc","--min_percent_cells",dest="min_percent_cells", type=float, default=None,  help="The minimum percent cells")
    parser.add_argument("-mcc","--min_count_cells",dest="min_count_cells", type=int, default=None,   help="The minimum count cells")
    parser.add_argument("-mnc","--min_num_cells",dest="min_num_cells", type=int, default=None,  help="The minimum number of cells")
    parser.add_argument("-ec","--expression_cutoff",dest="expression_cutoff", type=float, default=None, help="The expression cutoff")        
    parser.add_argument("-of","--of",dest="output_filename_prefix", default="StreamOutput",  help="output file name prefix")
                      




    args = parser.parse_args()
    
    input_filename = args.input_filename
    cell_label_filename = args.cell_label_filename
    cell_label_color_filename = args.cell_label_color_filename
    flag_norm = args.flag_norm
    flag_log2 = args.flag_log2
    output_folder = args.output_folder #work directory
    flag_remove_mt_genes = args.flag_remove_mt_genes
    min_count_genes = args.min_count_genes
    min_percent_cells = args.min_percent_cells
    min_percent_genes = args.min_percent_genes

    min_count_cells = args.min_count_cells
    min_num_cells = args.min_num_cells
    expression_cutoff = args.expression_cutoff
    output_filename_prefix = args.output_filename_prefix

    
    print('Starting mapping procedure...')
    if(output_folder==None):
        workdir_ref = os.path.join(os.getcwd(),'stream_result')
    else:
        workdir_ref = output_folder
    workdir = "./"

    if (input_filename.endswith('pkl')):
        adata = st.read(file_name=input_filename, file_format='pkl', workdir=workdir)
    else:
        adata = st.read(file_name=input_filename,workdir=workdir)
    print('Input: '+ str(adata.obs.shape[0]) + ' cells, ' + str(adata.var.shape[0]) + ' genes')

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if(cell_label_filename !=None):
        st.add_cell_labels(adata,file_name=cell_label_filename)
    else:
        st.add_cell_labels(adata)
    if(cell_label_color_filename !=None):   
        st.add_cell_colors(adata,file_name=cell_label_color_filename)
    else:
        st.add_cell_colors(adata)
        
    if (flag_norm):
        st.normalize_per_cell(adata)
 
    if (flag_log2):
        st.log_transform(adata, base=2)

    if (flag_remove_mt_genes):
        st.remove_mt_genes(adata)


    st.filter_cells(adata,min_pct_genes=min_percent_genes,min_count=min_count_genes, expr_cutoff=expression_cutoff)
    st.filter_genes(adata,min_num_cells=min_num_cells, min_pct_cells=min_percent_cells, min_count=min_count_genes, expr_cutoff=expression_cutoff)

    st.write(adata,file_name=(output_filename_prefix + '_stream_result.pkl'),file_path='./',file_format='pkl') 
    print('Output: '+ str(adata.obs.shape[0]) + ' cells, ' + str(adata.var.shape[0]) + ' genes')

    print('Finished computation.')

if __name__ == "__main__":
    main()
