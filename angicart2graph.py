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
    names = data[' name']
    names=  [str(i).strip().strip("\n") for i in names] 
    data[' name'] = names
    
    indeces = list(range(len(data[' name'])))
    name_index = {data[' name'][i]:indeces[i] for i in range(len(indeces))}
    n = len(indeces)
    adjmat = make_adjmat(data,name_index,indeces)

    features = np.array(pd.concat([data[' len(mm)'], data[' <r>_obs(mm)']], axis=1, keys=['length', 'radius']))
    
    features_dataf = pd.concat([data[' name'],data[' len(mm)'], data[' <r>_obs(mm)']], axis=1, keys=['id','length', 'radius'])

    ''' 
    graph = nx.from_numpy_array(adjmat)
    largest_cc = max(nx.connected_components(graph), key=len)
    graph = graph.subgraph(largest_cc).copy()
    adjmat = nx.to_numpy_array(graph)
    cc = np.array(list(largest_cc))
    features = np.array(features)[cc]
    '''
    return adjmat,features, name_index, features_dataf

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
    radius = tumor_radius*2
    
    near_tumor = []
    tumor_center = find_tumor_center(tumor)
    for key in progressbar.progressbar(name_coords):
        for point in name_coords[key]:
            if inSphere(point,tumor_center, radius):
                near_tumor.append(key)
                break
    return near_tumor


def get_vessels_in_tumor(tumor, name_coords):

    sensitivity = 0.5 #parameter for how much of a vessel needs to be in a tumor to consider bad.  = 0.5 means throw vessel if its half in the tumor.

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
    
    
    
    from mpl_toolkits.mplot3d import axes3d, Axes3D
    fig = plt.figure()
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
###############################
if __name__ == '__main__':
    #print(filenames)
    
    adjmat,features1, name_index, features1_dataf = process_angicart_output(filenames[1]) # return adjacency matrix and matrix of radius and length. same order as adjmat
    
    

    index_name = {}
    for key in name_index.keys():
        index_name[name_index[key]] = key
    
    in_tumor, near_tumor = process_vessel_tumor_coordinates(filenames[2],filenames[3], name_index) # returns names of vessels inside and around the tumor
    
    diseased = Diff(near_tumor, in_tumor)


    features1_dataf['labels'] = -1
    dis_inx= []
    for i in diseased:
        w = np.where(features1_dataf['id'] == i)[0][0]
        dis_inx.append(w)
    features1_dataf['labels'][dis_inx] = 1

    healthy = []
    while(len(healthy)<len(dis_inx)):
        r = np.random.randint(len(features1_dataf))
        if r not in dis_inx:
            healthy.append(r)
    features1_dataf['labels'][healthy] = 0

    #print(len(diseased))
    
    '''
    a = []
    b = []
    for i in near_tumor:
        a.append(name_index[i])
    for key in name_index.keys():
        if key not in in_tumor:
            b.append(name_index[key])

    
    viz_graph(adjmat,a,b)
    '''
    #print(filenames[0])

    features2 = pd.read_csv(join(input_path,filenames[0]), header = 0)
    names = features2.id
    names = [str(i).strip() for i in names]
    features2["id"] = names
    
    feats1_ids = features1_dataf["id"]
    feats2_ids = features2["id"]

    features_intersection = list(set(feats1_ids).intersection(set(feats2_ids)))
    
    feat1_inx= []
    for i in features_intersection:
        w = np.where(features1_dataf['id'] == i)[0][0]
        feat1_inx.append(w)
    
    feat2_inx= []
    for i in features_intersection:
        w = np.where(features2['id'] == i)[0][0]
        feat2_inx.append(w)

    #print(features1_dataf.iloc[feat1_inx].head())
    #features = pd.DataFrame([features1_dataf.iloc[feat1_inx],features2.iloc[feat2_inx]])
    
    feats1 = features1_dataf.iloc[feat1_inx].sort_values('id')
    feats2 = features2.iloc[feat2_inx].sort_values('id').drop("Unnamed: 0", axis = 1)
    #print(feats1.head())
    #print(feats2.head())
    #features2= features2.sort_values('indx')

    features = pd.merge(feats1,feats2,on='id')

    

    

    name = list(name_index.keys())
    diff = Diff(name,features_intersection)

    diff_inds = [name_index[i] for i in diff]

    graph = nx.from_numpy_array(adjmat)
    for n in diff_inds:
        graph.remove_node(n)
    adjmat = nx.to_numpy_array(graph)

    features.to_csv(join(input_path,folder_name+"_feature_matrix.csv"))
    adj_df = pd.DataFrame.from_records(adjmat)
    adj_df.to_csv(join(input_path,folder_name+"_adj_matrix.csv"))



    '''
    #print(features2.head())
    included= [str(i).strip().strip("\n") for i in list(features2['id'])] 
    names=  [str(i).strip().strip("\n") for i in list(name_index.keys())] 
    diseased=  [str(i).strip().strip("\n") for i in diseased]
    d = Diff(names,included)
    features_intersection = list(set(names).intersection(set(included)))
    print(d)
    name_index_clean ={}
    count = 0
    for key in name_index.keys():
        name = names[count]
        count+=1
        
        name_index_clean[name] = name_index[key]

    d = Diff(names,features_intersection)
    to_remove = [i for i in d] #list(set(list(in_tumor) + list(d)))
    to_remove =  list(set([str(i).strip().strip("\n") for i in to_remove]))
    adjmat,features1 = remove_vessels_by_name(adjmat,features1,name_index_clean,to_remove) # removes vessels from graph and feature matrix (used for vessels in tumor because invalid)
    features2[features2.columns[1]]=included

    diff = Diff(included,features_intersection)
    dr = []
   # w = np.where(features2.id in intersection)
    for i in diff:
        w = np.where(features2.id == i)[0][0]
        dr.append(w)
    dis_inds = []
    for i in diseased:
        w = np.where(features2.id == i)[0]
        if len(w)>0:
            dis_inds.append(w[0])

    print(len(features2))
    features2 = features2.drop(dr)
    print(len(features2))
    print(len(features1))
    inds = []
    for i in features2.id:
        if(i in name_index_clean):
            inds.append(name_index_clean[i])

    #print(features2.head())
    features2["len"] = features1[:,0] 
    features2["rad"] = features1[:,1]  
    features2["indx"] = inds
    features2 = features2.drop(["Unnamed: 0", "id"], axis = 1)
    features2= features2.sort_values('indx')
    features2["indx"] = list(range(len(features2["indx"])))
    features2["labels"] = -1
    features2["labels"].loc[dis_inds] = 1
    
   # print(features2["labels"][dis_inds])

    inds = list(range(len(features2)))
    healthy = []
    count = 0
    while(len(healthy)<len(dis_inds)):
        r = np.random.randint(len(inds))
        if r not in dis_inds:
            healthy.append(r)
            
    print(healthy)
    
    features2["labels"].loc[healthy] = 0
    #features2 = np.array(features2)
    #print(np.where(features2.labels == 1))
    features2.to_csv(join(input_path,folder_name+"_feature_matrix.csv"))
    adj_df = pd.DataFrame.from_records(adjmat)
    adj_df.to_csv(join(input_path,folder_name+"_adj_matrix.csv"))
    '''