import pandas as pd
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import FunctionTransformer
from os import listdir
from os.path import isfile, join, split
from os import walk

def load_all_data():
    #input_path = "D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\lung_adjmatrices_featurematrices\\"
    #input_path = "D:\work\Classes\Tufts\MLonGraphs\Project\data\lung_adjmatrices_featurematrices"

    input_path = "data/lung_adjmatrices_featurematrices/"

    to_drop = ["Unnamed: 0","id","labels"]
    max_graph_size = 0
    feat_size = 0

    f = []
    for (dirpath, dirnames, filenames) in walk(input_path):
        for filename in filenames:
            if '.tsv' in filename or '.csv' in filename:
                f.append(filename)
        break

    filenames = f
    filenames.sort()

    data = []
    
    for i in range(0,len(filenames)-1,2):
        adj = pd.read_csv(join(input_path,filenames[i]))
        features = pd.read_csv(join(input_path,filenames[i+1]))


        features["labels"][features["labels"] == -1] = -1 # set to -1 for semi supervised and 0 for naive supervised
        y = np.array(features["labels"])

        #A = lil_matrix(np.array(adj, dtype = "int32")[:,1:])
        features = features.drop(to_drop, axis = 1)
        
        pipe = Pipeline([
            
            ('scale', StandardScaler()),
            
        ])


        features = pipe.fit_transform(features)
        A = np.array(adj, dtype = "int32")[:,1:]
        X = np.array(features)
        data.append((A,X,y))
        feat_size = X.shape[-1]
        if A.shape[0]> max_graph_size:
            max_graph_size = A.shape[0]

    fixed_data = []
    for datum in data:
        A = datum[0]
        X = datum[1]
        y= datum[2]

        Abig = np.zeros((max_graph_size,max_graph_size))
        Abig = np.identity(max_graph_size)
        Abig[0:A.shape[0],0:A.shape[1]] = A

        Xbig = np.zeros((max_graph_size,feat_size))
        Xbig[0:X.shape[0],0:X.shape[1]] = X

        ybig = -1*np.ones((max_graph_size,))
        ybig[0:y.shape[0]] = y
        fixed_data.append((Abig,Xbig,ybig))
        
    return fixed_data, max_graph_size,feat_size


def load_data():
    input_path = "D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\lung_adjmatrices_featurematrices\\"
    #input_path = "data/lung_adjmatrices_featurematrices/"

    to_drop = ["Unnamed: 0","id","labels"]


    adj = pd.read_csv(join(input_path,"lung01_415_adj_matrix.csv"))
    features = pd.read_csv(join(input_path,"lung01_415_feature_matrix.csv"))

    features["labels"][features["labels"] == -1] = -1 # set to -1 for semi supervised and 0 for naive supervised
    y = np.array(features["labels"])

    #A = lil_matrix(np.array(adj, dtype = "int32")[:,1:])
    features = features.drop(to_drop, axis = 1)
    
    pipe = Pipeline([
        
        ('scale', StandardScaler()),
        
    ])


    features = pipe.fit_transform(features)
    A = np.array(adj, dtype = "int32")[:,1:]
    X = np.array(features)
    
    return A,X,y, A.shape[0], X.shape[-1]
