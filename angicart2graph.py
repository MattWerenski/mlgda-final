import argparse
from os import listdir
from os.path import isfile, join, split
from os import walk
import pandas as pd
import numpy as np

import networkx as nx
import matplotlib.pyplot as plt
import progressbar
from scipy.spatial import Delaunay
from scipy.spatial.distance import pdist, squareform




def make_adjmat(data, name_index, indeces):

    n = len(indeces)
    adjmat = np.zeros((n,n))

    for i in range(n):
        row = data.iloc[i]
        curr = name_index[row.iloc[0]] 
        children = []
        
        child_i = 7
        child = str(row.iloc[child_i]).strip()
        while not child == 'nan':
            
            children.append(name_index[child])
            child_i +=1
            if child_i >= len(row):
                break
            child = str(row.iloc[child_i]).strip()

        for child in children:
            adjmat[curr][child] = 1    
            adjmat[child][curr] = 1   
    return adjmat

def process_angicart_output(filename):
    
    data = pd.read_csv(join(input_path,filename), sep = '\t')
    data = data.sort_values(" name")
    
    names = data[' name']
    names=  [str(i).strip().strip("\n") for i in names] 
    data[' name'] = names
    to_drop = ['.04','.05','.06','.07','.08','.09']
    inds = []
    ind = 0
    for name in names:
        for d in to_drop:
            if d in name:
                inds.append(ind)
        ind+=1
    print(len(data))
    data = data.drop([data.index[inds[i]] for i in range(len(inds))])
    data = data.reindex(range(len(data)))
    print(len(data))
    print(data.head())
    indeces = list(range(len(data[' name'])))
    name_index = {data[' name'][i]:indeces[i] for i in range(len(indeces))}
    n = len(indeces)
    adjmat = make_adjmat(data,name_index,indeces)

    adjmat_df = pd.DataFrame(adjmat)
    adjmat_df["id"] = data[' name']
    
    adjmat_df_ = adjmat_df.drop("id", axis = 1)
    adjmat = adjmat_df_.to_numpy()
    

    features = np.array(pd.concat([data[' len(mm)'], data[' <r>_obs(mm)']], axis=1, keys=['length', 'radius']))
    
    features_dataf = pd.concat([data[' name'],data[' vol(cu.mm)'],data[' len(mm)'], data[' <r>_obs(mm)']], axis=1, keys=['id','vol','length', 'radius'])


    ''' 
    graph = nx.from_numpy_array(adjmat)
    largest_cc = max(nx.connected_components(graph), key=len)
    graph = graph.subgraph(largest_cc).copy()
    adjmat = nx.to_numpy_array(graph)
    cc = np.array(list(largest_cc))
    features = np.array(features)[cc]
    '''
    return adjmat,features, name_index, features_dataf, adjmat_df

def get_tumor_points(tumor_data):
    tumor = []
    for i in range(len(tumor_data)):
        point = np.array((tumor_data.iloc[i][0],tumor_data.iloc[i][1],tumor_data.iloc[i][2] +25))
        tumor.append(point)
    tumor = np.array(tumor)
    return tumor

def get_vessel_points(vessel_data, name_index):
    name_coords = {k: [] for k in name_index.keys()}
    for i in progressbar.progressbar(range(len(vessel_data))):
        point = np.array((vessel_data.iloc[i]['x'],vessel_data.iloc[i]['y'],vessel_data.iloc[i]['z']))
        name =  str(vessel_data.iloc[i]['nodeid']).strip()
        
        if name in name_coords.keys():
            #print(name)
            name_coords[name].append(point)
            
    for key in progressbar.progressbar(name_coords.keys()):
        name_coords[key] = np.array(name_coords[key])
        
    return name_coords

def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0


def inSphere(point, ref, radius):

    diff = np.subtract(point, ref)

    dist = np.sum(np.power(diff, 2))

    return dist < radius ** 2

def get_vessels_near_tumor(tumor, name_coords):
    D = pdist(tumor)
    D = squareform(D)
    tumor_radius = np.nanmax(D)/2
    radius = tumor_radius*1.2
    
    near_tumor = []
    tumor_center = find_tumor_center(tumor)
    for key in progressbar.progressbar(name_coords):
        for point in name_coords[key]:
            if inSphere(point,tumor_center, radius):
                near_tumor.append(key)
                break
    return near_tumor


def get_vessels_in_tumor(tumor, name_coords):

    sensitivity = 1.1 #parameter for how much of a vessel needs to be in a tumor to consider bad.  = 0.5 means throw vessel if its half in the tumor.

    tumor_hull = Delaunay(tumor)
    in_tumor = []
    for key in progressbar.progressbar(name_coords):
        if len(name_coords[key])>0:
            a = in_hull(name_coords[key],tumor_hull)
            if sum(a) >= len(a)*sensitivity:
                in_tumor.append(key)
    return in_tumor

def find_tumor_center(tumor):
    center = np.array([np.mean(tumor.T[0]),np.mean(tumor.T[1]),np.mean(tumor.T[2])])
    
    return center
def process_vessel_tumor_coordinates(vessel_file,tumor_file, name_index):
    vessel_data = pd.read_csv(join(input_path,vessel_file), sep = ',')
    tumor_data = pd.read_csv(join(input_path,tumor_file), sep = ',')

    tumor = get_tumor_points(tumor_data)
    name_coords = get_vessel_points(vessel_data,name_index)

    in_tumor = get_vessels_in_tumor(tumor,name_coords) #to be thrown away because false reconstructions inside solid mass
    near_tumor = get_vessels_near_tumor(tumor,name_coords) # to be labeled diseased (need to remove in_tumor first)
    
    '''
    #print(near_tumor)
    from mpl_toolkits.mplot3d import axes3d, Axes3D
    fig = plt.figure(1)
    ax = Axes3D(fig)
    ax.plot3D(tumor.T[0], tumor.T[1], tumor.T[2], 'gray')
    
    for key in name_coords.keys():
        ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'blue')
    
    for t in near_tumor:
        ax.plot3D(name_coords[t].T[0],name_coords[t].T[1],name_coords[t].T[2], color = 'green')
    for t in in_tumor:
        ax.plot3D(name_coords[t].T[0],name_coords[t].T[1],name_coords[t].T[2], color = 'red')
   

    fig = plt.figure(2)
    ax = Axes3D(fig)
    ax.plot3D(tumor.T[0], tumor.T[1], tumor.T[2], 'gray')
    

    #for key in name_coords.keys():
    #    ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'blue')
    
    for t in near_tumor:
        ax.plot3D(name_coords[t].T[0],name_coords[t].T[1],name_coords[t].T[2], color = 'green')
    for t in in_tumor:
        ax.plot3D(name_coords[t].T[0],name_coords[t].T[1],name_coords[t].T[2], color = 'red')
    #ax.plot3D([tumor_center[0],tumor_center[0]+0.5],[tumor_center[1],tumor_center[1]+0.5],[tumor_center[2],tumor_center[2]+0.5], color = 'green')

    plt.show()
    '''
    return in_tumor, near_tumor
    
    
def remove_vessels_by_name(adjmat,features,name_index,rem):
    graph = nx.from_numpy_array(adjmat)
    rem_nodes = []
    for name in rem:
        rem_nodes.append(name_index[name])
    
    for node in rem_nodes:
        graph.remove_node(node)
    new_features = []
    for i in range(len(features)):
        if i not in rem_nodes:
            new_features.append(features[i])

    new_adjmat = nx.to_numpy_array(graph)
    new_features = np.array(new_features)

    return new_adjmat,new_features


def viz_graph(adjmat, nodelist1,nodelist2):
    
    a= sum(adjmat)

    graph = nx.from_numpy_array(adjmat)
    #largest_cc = max(nx.connected_components(graph), key=len)
    #graph = graph.subgraph(largest_cc).copy()
    #A_lcc = nx.to_numpy_array(graph)
    #a= sum(A_lcc)
    #plt.hist(a)
    
    pos = nx.kamada_kawai_layout(graph)
    nx.draw_networkx_nodes(graph,pos, node_size=10, node_color = 'b', nodelist = nodelist2)
    nx.draw_networkx_nodes(graph,pos, node_size=10, node_color = 'r', nodelist = nodelist1)
    nx.draw_networkx_edges(graph,pos)
    plt.show()

def Diff(li1, li2): 
    return (list(set(li1) - set(li2))) 

from sklearn.neighbors import NearestNeighbors
def add_spatial_connections(vessel_file,adjmat_ids, valid_ids= None, adjmat = None): # This makes spacial edges

    if valid_ids is None: #Ignore please, go to else
        vessel_data = pd.read_csv(join(input_path,vessel_file), sep = ',')
        ids = np.array(adjmat_ids["id"])
        name_index = {ids[i]:i for i in range(len(ids))}
        
        name_points = get_vessel_points(vessel_data,name_index)
        
        coms = dict()
        
        for name in name_points.keys():
            if name in ids:
                vessel = name_points[name]
                if len(vessel) !=0:
                
                    com = np.mean(vessel, axis = 0)
                    coms[name] = com
                
        coms_vals = np.array(list(coms.values()))
        K = 4
        neigh = NearestNeighbors(K, 1, algorithm = 'brute',metric = "mahalanobis", metric_params= {'V':np.cov(coms_vals)})
        print(coms_vals.shape)
        neigh.fit(coms_vals)
        edges = []
        for name in progressbar.progressbar(coms.keys()):
            neighbs = neigh.kneighbors([coms[name]], K, return_distance=False)
            for i in neighbs[0][1:]:
                #print(name)
                edges.append((name,ids[i]))
                row = adjmat_ids[adjmat_ids['id'] == name]
                #print(np.where(np.array(row) == 1))
                
                
                adjmat_ids.loc[adjmat_ids['id'] == name,i] = 1
                l = adjmat_ids.index[adjmat_ids['id'] == name].tolist()
                
                adjmat_ids.loc[i,l] = 1
                #print(adjmat_ids[adjmat_ids['id'] == name][i])
    else:
        # This is the working correct thing
        vessel_data = pd.read_csv(join(input_path,vessel_file), sep = ',')
        ids = np.array(vessel_data["nodeid"])
        name_index = {str(ids[i]):i for i in range(len(ids))}
        print(type(valid_ids[0]))
        name_points = get_vessel_points(vessel_data,name_index)
        
        coms = dict()
        
        for name in valid_ids:
            vessel = name_points[name]
            if len(vessel) !=0:
            
                com = np.mean(vessel, axis = 0)
                coms[name] = com
                
        coms_vals = np.array(list(coms.values()))
        neigh = NearestNeighbors(4, 1) # 4 neighbors because first is self. I left it in there because why not. 
        print(coms_vals.shape)
        neigh.fit(coms_vals)

        edges = []
        for name in progressbar.progressbar(coms.keys()):
            neighbs = neigh.kneighbors([coms[name]], 4, return_distance=False)
            for i in neighbs[0][1:]:
                #print(name)
                edges.append((name,valid_ids[i]))
        return edges #used in make_valid_adj_matrix





    return adjmat_ids



def make_valid_adj_matrix(filename, valid_ids):
    data = pd.read_csv(join(input_path,filename), sep = '\t')
    data = data.sort_values(" name")
    
    names = data[' name']
    names=  [str(i).strip().strip("\n") for i in names] 
    data[' name'] = names

    n = len(data)
    adjlist = dict()

    for i in range(n):
        row = data.iloc[i]
        curr = row.iloc[0]
        children = []
        
        child_i = 7
        child = str(row.iloc[child_i]).strip()
        while not child == 'nan':
            
            if child in valid_ids:
                children.append(child)
            child_i +=1
            if child_i >= len(row):
                break
            child = str(row.iloc[child_i]).strip()

        if curr in valid_ids:
            adjlist[curr]  = children
    

    n = len(list(adjlist.keys()))
    indeces = list(range(n))
    keys = list(adjlist.keys())
    name_index = {keys[i]:indeces[i] for i in range(n)}

    vals = list(adjlist.values())
    adjmat = np.zeros((n,n))
    for i in range(n):
        for c in vals[i]:
            j = name_index[c]
            adjmat[i][j]=1
            adjmat[j][i]=1


    edges = add_spatial_connections(filenames[2],None,valid_ids=valid_ids)
    for v,u in edges:
        i = name_index[v]
        j = name_index[u]
        adjmat[i][j]=1
        adjmat[j][i]=1


    
    return adjmat





###############################
if __name__ == '__main__':
    #print(filenames)
    parser = argparse.ArgumentParser(description='Covert Angicart++ output to adjacency matrix and file of features.')
    parser.add_argument('path_to_files', metavar='P', type=str, nargs=1,
                        help='A path to the folder containing angicart output files.')

    args = parser.parse_args()

    input_path = args.path_to_files[0]

    folder_name = split(split(input_path)[0])[1]
    print(folder_name)

    f = []
    for (dirpath, dirnames, filenames) in walk(input_path):
        for filename in filenames:
            if '.tsv' in filename or '.csv' in filename:
                f.append(filename)
        break

    filenames = f
    
    adjmat,features1, name_index, features1_dataf, adjmat_ids = process_angicart_output(filenames[1]) # return adjacency matrix and matrix of radius and length. same order as adjmat
    adjmat_ids = add_spatial_connections(filenames[2],adjmat_ids) # ignore
    

    index_name = {}
    for key in name_index.keys():
        index_name[name_index[key]] = key
    
    in_tumor, near_tumor = process_vessel_tumor_coordinates(filenames[2],filenames[3], name_index) # returns names of vessels inside and around the tumor
    
    diseased = Diff(near_tumor, in_tumor)

    subtree = near_tumor[0].split('.')[1]
    print("-----------------------------")
    print(subtree)
    features1_dataf['labels'] = -1
    dis_inx= []
    labels = -1*np.ones(len(features1_dataf))
    ids = np.array(features1_dataf['id'])
    for i in range(len(ids)):
        if ids[i] in diseased:
            labels[i] = 1
    
    for i in diseased:
        w = np.where(features1_dataf['id'] == i)[0][0]
        
        #w = features1_dataf['id'][features1_dataf['id'] == i]
        dis_inx.append(w)
    features1_dataf['labels']=labels

    healthy = []
    healthy_ids = []
    while(len(healthy)<len(diseased)*3):
        r = np.random.randint(len(features1_dataf))
        if features1_dataf['id'][r] not in diseased:
            if subtree not in str(features1_dataf['id'][r]):
                healthy.append(r)
                healthy_ids.append(features1_dataf['id'][r])
    features1_dataf['labels'][healthy] = 0
    print(near_tumor)
    print(healthy_ids)

    #print(filenames[0])

    features2 = pd.read_csv(join(input_path,filenames[0]), header = 0)
    names = features2.id
    names = [str(i).strip() for i in names]
    features2["id"] = names
    
    feats1_ids = features1_dataf["id"]
    feats2_ids = features2["id"]

    features_intersection = list(set(feats1_ids).intersection(set(feats2_ids)))

    #feats1= features1_dataf[features1_dataf['id'] == features_intersection]
    #feats2= features2[features2['id'] == features_intersection]

    features = pd.merge(features1_dataf,features2,on='id').drop("Unnamed: 0", axis = 1)
    
    valid_ids = np.array(features['id'])

    
    print(features.head())
   
    adjmat = make_valid_adj_matrix(filenames[1],valid_ids)
    adjmat = pd.DataFrame.from_records(adjmat)
    #adjmat.to_csv("test.csv")
    features.to_csv(join("D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\lung_adjmatrices_featurematrices",folder_name+"_feature_matrix.csv"))
    adj_df = pd.DataFrame.from_records(adjmat)
    adjmat.to_csv(join("D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\lung_adjmatrices_featurematrices",folder_name+"_adj_matrix.csv"))


