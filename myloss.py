import tensorflow as tf
from tensorflow import keras
import spektral
import keras
from spektral.layers import GraphSageConv, GraphConv
from keras.layers import Input, Dense
from keras.models import Model
from keras.optimizers import Adam
from keras.layers import Layer
from keras.activations import relu
import numpy as np
import keras.backend as K

def semi_supervised(y,preds):
    where = keras.backend.cast(y > -1 , dtype = 'int32')
    y_pred = K.gather(preds,where)
    y_true =  K.gather(y,where)
    loss = K.binary_crossentropy(y_true, y_pred)
    return loss

def semi_supervised_acc(y,preds):
    where = keras.backend.cast(y > -1 , dtype = 'int32')
    y_pred = K.gather(preds,where)
    y_true =  K.gather(y,where)
    acc = keras.metrics.binary_accuracy(y_true, y_pred)
    #acc = keras.backend.cast(y_true == y_pred , dtype = 'int32')
    #acc = K.mean(acc)
    return acc