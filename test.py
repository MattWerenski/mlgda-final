import tensorflow as tf
from tensorflow import keras 
import pandas as pd
import numpy as np
from model import SimpleModel, DenseModel
from scipy.sparse import lil_matrix
from load_data import *



from tensorflow import set_random_seed
set_random_seed(2)
np.random.seed(2)

#data, max_graph_size, feat_size = load_all_data()
data, max_graph_size, feat_size = load_all_data_two_matrices()
#A,X,y = data[0]


#A,X,y, max_graph_size, feat_size = load_data()



adjsv = []
adjss = []
feats = []
labels = []

for i in [0,2,3,4,5]:
    adjsv.append(data[i][0])
    adjss.append(data[i][1])
    feats.append(data[i][2])
    labels.append(data[i][3])

print('-------- SHAPES ----------')
print('max graph ',max_graph_size)
print(np.concatenate(labels).shape)
print(np.concatenate(adjsv).shape)
print(np.concatenate(feats).shape)



#print(max_graph_size, feat_size)#
model = SimpleModel(max_graph_size, feat_size)
#model = DenseModel(max_graph_size, feat_size)

Av = np.concatenate(adjsv)
As = np.concatenate(adjss)
Xs = np.concatenate(feats)
ys = np.concatenate(labels)


model.fit(As,Xs,ys)
#model.fit_many(data[1:])

A,X,y = data[1]
score1, conf, auc, f1 = model.evaluate(A,X,y)


print(score1, conf, auc, f1)
preds= np.round(model.predict(A,X))


input_path = "D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\lung_adjmatrices_featurematrices\\"
features = pd.read_csv(join(input_path,"lung01_041_feature_matrix.csv"))

df= pd.DataFrame.from_records(preds)
df['id'] = features['id']
df.to_csv("predictions_lung01_041.csv")
#print(score2, conf2)
