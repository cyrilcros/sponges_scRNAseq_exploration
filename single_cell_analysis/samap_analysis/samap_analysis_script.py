#################
### run SAMap ###
#################
from samalg.utilities import convert_annotations
from samap.analysis import get_mapping_scores
from samap import q
import rpy2.robjects.pandas2ri
import rpy2.robjects.numpy2ri
from rpy2.robjects import r
from samap.mapping import SAMAP
from samap.analysis import get_mapping_scores, GenePairFinder
from samalg import SAM
import samalg.utilities as ut
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def plot_sankey(CSIM,A1,B1,nid,dom, w=1000,h=1500,remove_zeros=True,thr=1.0,CSIMthr=None,group=None,posy=None,max_dummies=100,ADD_REF_LINK=False,make_pdf=False):
    if CSIMthr is None:
        CSIMthr = CSIM.copy()
    X = CSIM.copy()
    X[CSIMthr<thr]=0
    Y = np.vstack((np.hstack((np.zeros((A1.size,A1.size)),X)),np.hstack((X.T*0,np.zeros((B1.size,B1.size))))))
    Aorig = A1.copy()
    Borig = B1.copy()
    if remove_zeros:
        idx = np.where(Y.sum(1)>0)[0]
        #idx = idx#idx[idx>=A1.size]
        Y2 = Y[idx,:]
        At=A1
        A1=A1[idx]
        idx = np.where(Y2.sum(0)>0)[0]
        Y2 = Y2[:,idx]
        B1=np.append(At,B1)[idx]
        Y2 = np.hstack((np.zeros((Y2.shape[0],Y2.shape[0])),Y2))

    As,Bs=Aorig.size,Borig.size
    Z = np.zeros((max(Aorig.size,Borig.size)*2,)*2)
    Z[:Y.shape[0],:Y.shape[1]]=Y
    if Aorig.size>Borig.size:
        FIRST=False
        dummy='dummy'+np.arange(min(max_dummies,Aorig.size-Borig.size)).astype('<U50').astype('object')
        Borig = np.append(Borig,dummy)
    elif Borig.size > Aorig.size:
        FIRST=True
        dummy = 'dummy'+np.arange(min(max_dummies,Borig.size-Aorig.size)).astype('<U50').astype('object')
        Aorig = np.append(Aorig,dummy)
    else:
        FIRST=False
    Y = Z
    x,y = Y.nonzero()
    data = Y[x,y]

    if FIRST:
        x[x>=As] += dummy.size
        y[y>=As] += dummy.size

    if ADD_REF_LINK:
        Aorig = np.append(Aorig,'dummya')
        Borig = np.append(Borig,'dummyb')
        x[x>=Aorig.size-1]+=1
        y[y>=Aorig.size-1]+=1
        data = np.append(data,1)
        x=np.append(x,Aorig.size-1)
        y=np.append(y,Aorig.size+Borig.size-1)

    LINKS=pd.DataFrame(list(zip(x,y,data)),columns=['source','target','value'])#.to_csv('/media/storage/SAMmap/'+nid+'_sankey_links.csv')
    #clu = np.append('A_'+clu1.astype('object'),'B_'+clu2.astype('object'))
    clu = np.append(Aorig,Borig)#[idx]
    NODES=pd.DataFrame(data = np.vstack((np.arange(clu.size),clu)).T,columns=['node','name'])#.to_csv('/media/storage/SAMmap/'+nid+'_sankey_nodes.csv')
    if group is not None:
        #'/media/storage/SAMmap/ZFXF_groups.csv'
        if type(group) is str:
            C = pd.read_csv(group,index_col=0)
        else:
            C=group
        clu2=clu.copy()
        filt = np.in1d(clu2,ut.search_string(clu,'dummy')[0],invert=True)
        clu2 = clu2[filt]
        vals = np.zeros(clu.size)-1
        vals[filt] = C.T[clu2].T.values.flatten()
        NODES['group'] = vals+1
    NODES['xPos'] = np.array([0]*Aorig.size+[1]*Borig.size)#[idx]
    if posy is not None:
        #'/media/storage/SAMmap/ZFXF_posy.csv'
        if type(posy) is str:
            C = pd.read_csv(posy, index_col=0)
        else:
            C = posy

        clu2=clu.copy()
        filt = np.in1d(clu2,ut.search_string(clu,'dummy')[0],invert=True)
        clu2 = clu2[filt]
        vals = np.zeros(clu.size)+10000
        vals[filt] = C.T[clu2].T.values.flatten()
        NODES['yPos'] = vals
        Anot = Aorig[np.in1d(Aorig,A1,invert=True)]
        Bnot = Borig[np.in1d(Borig,B1,invert=True)]
        NODES.iloc[np.where(np.in1d(NODES['name'].values,np.append(Anot,Bnot)))[0],-1] = NODES.iloc[np.where(np.in1d(NODES['name'].values,np.append(Anot,Bnot)))[0],-1].values.flatten()+1000
    sankey(NODES,LINKS,nid+'_sankey.html',dom,w,h)
    if make_pdf:
        pdfkit.from_file(nid+'_sankey.html', nid+'_sankey.pdf');
        
def compute_csim(sam3, key, X=None, n_top = 100):
    cl1 = q(sam3.adata.obs[key].values[sam3.adata.obs["batch"] == "batch1"])
    clu1,cluc1 = np.unique(cl1,return_counts=True)
    cl2 = q(sam3.adata.obs[key].values[sam3.adata.obs["batch"] == "batch2"])
    clu2,cluc2 = np.unique(cl2,return_counts=True)

    clu1s = q("batch1_" + clu1.astype("str").astype("object"))
    clu2s = q("batch2_" + clu2.astype("str").astype("object"))
    cl = q(
        sam3.adata.obs["batch"].values.astype("object")
        + "_"
        + sam3.adata.obs[key].values.astype("str").astype("object")
    )

    CSIM1 = np.zeros((clu1s.size, clu2s.size))
    if X is None:
        X = sam3.adata.obsp["connectivities"].copy()

    for i, c1 in enumerate(clu1s):
        for j, c2 in enumerate(clu2s):
            CSIM1[i, j] = np.append(
                np.sort(X[cl == c1, :][:, cl == c2].sum(1).A.flatten())[::-1][:n_top].mean(),
                np.sort(X[cl == c2, :][:, cl == c1].sum(1).A.flatten())[::-1][:n_top].mean(),
            ).max()
    CSIMth = CSIM1 / sam3.adata.uns['mdata']['knn_1v2'][0].data.size    
    s1 = CSIMth.sum(1).flatten()[:, None]
    s2 = CSIMth.sum(0).flatten()[None, :]
    s1[s1 == 0] = 1
    s2[s2 == 0] = 1
    CSIM1 = CSIMth / s1
    CSIM2 = CSIMth / s2
    CSIM = (CSIM1 * CSIM2) ** 0.5

    return CSIM, clu1, clu2, CSIMth

if __name__ == "__main__":
    id1 = 'sp'
    id2 = 'am'

    # passing in file names (SAMap will process the data with SAM and save the resulting objects to two `.h5ad` files.)
    fn1 = 'spongilla_run.h5ad' #processed data will be automatically saved to `/path/to/file/file1_pr.h5ad`
    fn2 = 'adult_amphimedon.h5ad' #processed data will be automatically saved to `/path/to/file/file2_pr.h5ad`
    sam1=SAM()
    sam2=SAM()
    sam1.load_data(fn1)
    sam2.load_data(fn2)

    sm = SAMAP(sam1,sam2,id1,id2,f_maps = 'maps/', save_processed=True)
    samap = sm.run(NH1=2,NH2=2,NUMITERS=6)


    ####################################
    ### plot gene expression overlap ###
    ####################################

    gene_pairs = ['sp_c100570-g1;am_Aqu2.1.22256_001',
    'sp_c100709-g1;am_Aqu2.1.39899_001',
    'sp_c100764-g2;am_Aqu2.1.21064_001',
    'sp_c100793-g2;am_Aqu2.1.34863_001',
    'sp_c100867-g1;am_Aqu2.1.34613_001',
    'sp_c100907-g1;am_Aqu2.1.35107_001',
    'sp_c101215-g1;am_Aqu2.1.31581_001',
    'sp_c101225-g1;am_Aqu2.1.39160_001',
    'sp_c101324-g1;am_Aqu2.1.33501_001',
    'sp_c101611-g1;am_Aqu2.1.41646_001',
    'sp_c101977-g1;am_Aqu2.1.34457_001',
    'sp_c102126-g1;am_Aqu2.1.39257_001',
    'sp_c102242-g1;am_Aqu2.1.37539_001',
    'sp_c102253-g1;am_Aqu2.1.36288_001',
    'sp_c102281-g3;am_Aqu2.1.18692_001',
    'sp_c102593-g1;am_Aqu2.1.21547_001',
    'sp_c102820-g1;am_Aqu2.1.31233_001',
    'sp_c103327-g3;am_Aqu2.1.40413_001',
    'sp_c103609-g2;am_Aqu2.1.43326_001',
    'sp_c103727-g1;am_Aqu2.1.22171_001',
    'sp_c103973-g1;am_Aqu2.1.43332_001',
    'sp_c104008-g2;am_Aqu2.1.42270_001',
    'sp_c104085-g1;am_Aqu2.1.33096_001',
    'sp_c104125-g2;am_Aqu2.1.33768_001',
    'sp_c104180-g1;am_Aqu2.1.40603_001',
    'sp_c104336-g2;am_Aqu2.1.32579_001',
    'sp_c104430-g1;am_Aqu2.1.40987_001',
    'sp_c104982-g1;am_Aqu2.1.36291_001',
    'sp_c110160-g1;am_Aqu2.1.44167_001',
    'sp_c110173-g1;am_Aqu2.1.34476_001',
    'sp_c110863-g1;am_Aqu2.1.43341_001',
    'sp_c2276-g1;am_Aqu2.1.40409_001',
    'sp_c39449-g1;am_Aqu2.1.37127_001',
    'sp_c40200-g1;am_Aqu2.1.28816_001',
    'sp_c77869-g1;am_Aqu2.1.34257_001',
    'sp_c78291-g1;am_Aqu2.1.18688_001',
    'sp_c78652-g1;am_Aqu2.1.39626_001',
    'sp_c79258-g1;am_Aqu2.1.43382_001',
    'sp_c79553-g1;am_Aqu2.1.41085_001',
    'sp_c79748-g1;am_Aqu2.1.23379_001',
    'sp_c82957-g1;am_Aqu2.1.18087_001',
    'sp_c85099-g1;am_Aqu2.1.33636_001',
    'sp_c85951-g1;am_Aqu2.1.43300_001',
    'sp_c87612-g1;am_Aqu2.1.31166_001',
    'sp_c90372-g1;am_Aqu2.1.41064_001',
    'sp_c92216-g1;am_Aqu2.1.32930_001',
    'sp_c93752-g1;am_Aqu2.1.32963_001',
    'sp_c94243-g1;am_Aqu2.1.35077_001',
    'sp_c94782-g1;am_Aqu2.1.38191_001',
    'sp_c95300-g1;am_Aqu2.1.36695_001',
    'sp_c95342-g1;am_Aqu2.1.22581_001',
    'sp_c95982-g1;am_Aqu2.1.38796_001',
    'sp_c96134-g1;am_Aqu2.1.27828_001',
    'sp_c97430-g1;am_Aqu2.1.38181_001',
    'sp_c97536-g1;am_Aqu2.1.41386_001',
    'sp_c97882-g1;am_Aqu2.1.27188_001',
    'sp_c97920-g1;am_Aqu2.1.27775_001',
    'sp_c98298-g1;am_Aqu2.1.34860_001',
    'sp_c98559-g1;am_Aqu2.1.36076_001',
    'sp_c98588-g1;am_Aqu2.1.14505_001',
    'sp_c98839-g2;am_Aqu2.1.13103_001',
    'sp_c99321-g1;am_Aqu2.1.32600_001',
    'sp_c99696-g1;am_Aqu2.1.30028_001',
    'sp_c99761-g1;am_Aqu2.1.28952_001',
    'sp_c99909-g1;am_Aqu2.1.33667_001']

    gene_pairs = [x.split(';') for x in gene_pairs]


    ut.create_folder('expression_umaps_pdfs')
    for I in range(len(gene_pairs)):
        print(I)
        a = gene_pairs[I][0]
        b = gene_pairs[I][1]
        ax = sm.plot_expression_overlap(a,b,s0=1,s1=3,s2=3,s3=10,thr=0.1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.figure.set_size_inches((6,6))
        ax.figure.savefig(f'expression_umaps_pdfs/{a};{b}.pdf',transparent=True,bbox_inches='tight')
        plt.close('all')

    ############################
    ### generate sankey plot ###
    ############################



    rpy2.robjects.pandas2ri.activate()
    rpy2.robjects.numpy2ri.activate()

    sankey = r('''
                 function(nodes,links,fname,dom,height,width){
                     #my_color <- 'd3.scaleOrdinal().domain([0,1,2,3]).range(["#24b324", "#2424b3","#b36b24","#b32424"])'
                     if( "yPos" %in% colnames(nodes)){
                         ypos <- "yPos"
                     }else{
                         ypos<-NULL
                     }
                     if( "group" %in% colnames(nodes)){
                         grp <- "group"
                     }else{
                         grp<-'name'
                     }
                         net<- sankeyD3::sankeyNetwork(Links = links, Nodes = nodes, NodePosX = 'xPos',
                             Source = 'source', highlightChildLinks=FALSE,dragY=TRUE,showNodeValues=FALSE,
                             Target = 'target', scaleNodeBreadthsByString=FALSE, align='none',
                             Value = 'value',
                             NodeID = 'name', NodePosY = ypos,
                             units = 'ew',#colourScale=my_color,
                             fontSize=8,fontFamily='Arial',
                             height=height,NodeGroup = grp,
                             width=width, xAxisDomain = dom)
                        sankeyD3::saveNetwork(net, fname,selfcontained=TRUE)
                 }
                 ''')
    r('''
    library(dplyr)
    library(sankeyD3)
    library(tidyr)
    ''')

    get_mapping_scores(sm,'CellTypeShort','celltype');
    A = pd.read_csv('cell_types_ordered.txt',sep='\t')
    group = convert_annotations(A['Family'].values.flatten())
    index = np.array(A['cell_type_shortname'])
    order = np.arange(index.size)
    n=sm.sam2.adata.obs['celltype'].cat.categories
    index = np.append(index,n)
    group = np.append(group,[group.max()+1]*n.size)
    order = np.append(order,[order.max()+1]*n.size)

    pd.DataFrame(data = group,index=index).to_csv('spam_group.csv')
    pd.DataFrame(data = order,index=index).to_csv('spam_order.csv')

    CSIMx,ax,bx,CSIMthrx = compute_csim(sm.samap,'CellTypeShort;celltype_mapping_scores', n_top = 1000000)
    nid = 'SPAM' 
    s3=sm.samap
    thr=0.15

    CSIMthrx[CSIMthrx<thr]=0

    w = 91 / 25.4 * 96 * 1.25**2.5
    h = 120 / 25.4 * 96 * 1.25
    CSIMx[CSIMx<0.]=0

    try:
        n = plot_sankey(CSIMthrx,q([x[3:] for x in ax]),q([x[3:] for x in bx]),nid,['Spongilla','Amphimedon'],ADD_REF_LINK=True,remove_zeros=False,thr=0,w=h,h=w,CSIMthr=CSIMthrx,make_pdf=True,max_dummies=1,
                    group='spam_group.csv',posy='spam_order.csv')
    except NameError:
        pass;

    #####################
    ### generate UMAP ###
    #####################

    ax = sm.scatter(c1='#000098',c2='#ffb900',s2=10,alpha2=0.8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.figure.set_size_inches((7,7))
    ax.figure.savefig('spam_umap.pdf',transparent=True,bbox_inches='tight')