import pandas as pd
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.preprocessing import FunctionTransformer
from os import listdir
from os.path import isfile, join, split
from os import walk
import os 

def load_all_data():
    #input_path = "D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\lung_adjmatrices_featurematrices\\"
    input_path = "D:\work\Classes\Tufts\MLonGraphs\Project\data\lung_adjmatrices_featurematrices"

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
            ('poly', PolynomialFeatures(degree =2))
            
        ])

        
        
        features["radius"]=features["radius"].replace({0:np.median(features.radius)})
        features["length"]=features["length"].replace({0:np.median(features.length)})
        features["TC"]=features["TC"].replace({0:np.median(features.TC)})
        features["AC"]=features["AC"].replace({0:np.median(features.AC)})
        features["TT"]=features["TT"].replace({0:np.median(features.TT)})
        features["MC"]=features["MC"].replace({0:np.median(features.MC)})
        features["MT"]=features["MT"].replace({0:np.median(features.MT)})
        features["TCC"]=features["TCC"].replace({0:np.median(features.TCC)})
        features["ACC"]=features["ACC"].replace({0:np.median(features.ACC)})
        features["TTsq"]=features["TTsq"].replace({0:np.median(features.TTsq)})
        features["dist_met"]=features["dist_met"].replace({0:np.median(features.dist_met)})
        features["inflec_met"]=features["inflec_met"].replace({0:np.median(features.inflec_met)})
        features["soam"]=features["soam"].replace({0:np.median(features.soam)})

        features["mean_ratio"]= np.abs(features["mean_ratio"])
       # features["mean_ratio"]= np.maximum(features["mean_ratio"],1000)

        features["mean_ratio"]=features["mean_ratio"].replace({0:np.median(features.max_ratio)})
        features["max_ratio"]= np.abs(features["max_ratio"])
        #features["max_ratio"]= np.maximum(features["max_ratio"],1000)
        features["max_ratio"]=features["max_ratio"].replace({0:np.median(features.max_ratio)})

        features['volume'] = features.radius * features.length
        features.length = np.log(features.length)
        features.radius = np.log(features.radius)
        features.volume = np.log(features.volume)
        features.TC = np.log(features.TC)
        features.AC = np.log(features.AC)
        features.TT = np.log(features.TT)
        features.MC = np.log(features.MC)
        features.MT = np.log(features.MT)
        features.TCC = np.log(features.TCC)
        features.ACC = np.log(features.ACC)
        features.TTsq = np.log(features.TTsq)
        features.dist_met = np.log(features.dist_met)
        features.inflec_met = np.log(features.inflec_met)
        features.soam = np.log(features.soam)
        features.mean_ratio = np.clip(features.mean_ratio,0,1)
        features.max_ratio= np.clip(features.max_ratio,0,1)
        features.mean_ratio = np.log(features.mean_ratio)
        features.max_ratio= np.log(features.max_ratio)

        features = features.drop(['mean_ratio','max_ratio'],axis = 1)

        features = features[['radius', 'length',"vol",'volume','dist_met','soam','inflec_met']]


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

        deg_mat_inv_half = np.zeros(A.shape)
        for j in range(len(A)):
            deg_mat_inv_half[j][j] = sum(A[j])**(-0.5)

        A = np.matmul(deg_mat_inv_half,A)
        A = np.matmul(A,deg_mat_inv_half) 

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
    features.length = np.log(features.length)
    features.radius = np.log(features.radius)
    features.TC = np.log(features.TC)
    features.AC = np.log(features.AC)
    features.TT = np.log(features.TT)
    features.MC = np.log(features.MC)
    features.MT = np.log(features.MT)
    features.TCC = np.log(features.TCC)
    features.ACC = np.log(features.ACC)
    features.TTsq = np.log(features.TTsq)
    features.dist_met = np.log(features.dist_met)
    features.inflec_met = np.log(features.inflec_met)
    features.soam = np.log(features.soam)

    pipe = Pipeline([
        
        ('scale', StandardScaler()),
        
    ])


    features = pipe.fit_transform(features)
    A = np.array(adj, dtype = "int32")[:,1:]
    X = np.array(features)
    
    return A,X,y, A.shape[0], X.shape[-1]