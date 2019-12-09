import tensorflow as tf
from tensorflow import keras 
import pandas as pd
import numpy as np
from model import SimpleModel
from scipy.sparse import lil_matrix
from load_data import *

#data, max_graph_size, feat_size = load_all_data()
#A,X,y = data[0]


A,X,y, max_graph_size, feat_size = load_data()

print('-------- SHAPES ----------')
print(A.shape, X.shape, y.shape)
print(max_graph_size, feat_size)
model = SimpleModel(max_graph_size, feat_size)

model.fit(A,X,y)

score = model.evaluate(A,X,y)

print(score)

