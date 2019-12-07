import pandas as pd
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import FunctionTransformer

from os.path import join

def load_data():
    #input_path = "D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\lung_adjmatrices_featurematrices\\"
    input_path = "data/lung_adjmatrices_featurematrices/"

    to_drop = ["Unnamed: 0","id","labels"]


    adj = pd.read_csv(join(input_path,"lung01_024_adj_matrix.csv"))
    features = pd.read_csv(join(input_path,"lung01_024_feature_matrix.csv"))

    features["labels"][features["labels"] == -1] = 0 # set to -1 for semi supervised and 0 for naive supervised
    y = np.array(features["labels"])

    #A = lil_matrix(np.array(adj, dtype = "int32")[:,1:])
    features = features.drop(to_drop, axis = 1)
    
    pipe = Pipeline([
        
        ('scale', StandardScaler()),
        
    ])


    features = pipe.fit_transform(features)
    A = np.array(adj, dtype = "int32")[:,1:]
    X = np.array(features)
    
    return A,X,y
