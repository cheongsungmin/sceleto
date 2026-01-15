from __future__ import annotations
from typing import Optional, Sequence, Union, List

def us(adata,gene,groups=None, show=False, exclude =None,figsize=None,**kwargs):
    """
    * 03/10/2022
    Create a umap using a list of genes.

    adata:AnnData,         REQUIRED | AnnData object.
    gene:list/str,         REQUIRED | List of genes to use for UMAP. A coma seperated string can be used instead of a list
    groups:str,        NOT REQUIRED | Restrict to a few categories in categorical observation annotation
    show:boolean,      NOT REQUIRED | Show the plot. Default = False.
    exclude:list,      NOT REQUIRED | List of genes to exclude. 
    figsize:float,     NOT REQUIRED | Figure size.
    """
    import scanpy as sc
    from matplotlib import rcParams

    if figsize:
        rcParams['figure.figsize'] = figsize
    if ',' in gene:
        gene = gene.split(',')
    if groups:
        sc.pl.umap(adata,color=gene,color_map='OrRd',groups=groups, show=show, **kwargs)
    else:
        if exclude: # list to not show
            groups = [x for x in set(adata.obs[gene]) if x not in exclude]
            sc.pl.umap(adata,color=gene,color_map='OrRd',groups=groups, show=show, **kwargs)
        else:
            sc.pl.umap(adata,color=gene,color_map='OrRd',show=show, **kwargs)
    rcParams['figure.figsize'] = [5,5]