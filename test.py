import tensorflow as tf
from tensorflow import keras 
import pandas as pd
import numpy as np
from model import SimpleModel
from scipy.sparse import lil_matrix
from load_data import *

data, max_graph_size, feat_size = load_all_data()
#A,X,y = data[0]


#A,X,y, max_graph_size, feat_size = load_data()

adjs = []
feats = []
labels = []

for i in range(4):
    adjs.append(data[i][0])
    feats.append(data[i][1])
    labels.append(data[i][2])

print('-------- SHAPES ----------')
print('max graph ',max_graph_size)
print(np.concatenate(labels).shape)
print(np.concatenate(adjs).shape)
print(np.concatenate(feats).shape)



#print(max_graph_size, feat_size)
model = SimpleModel(max_graph_size, feat_size)

As = np.concatenate(adjs)
Xs = np.concatenate(feats)
ys = np.concatenate(labels)


model.fit(As,Xs,ys)
#model.fit_many(data[1:])

A,X,y = data[0]
score = model.evaluate(A,X,y)

print(score)

