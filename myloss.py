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
    where = keras.backend.cast(y > -1 , dtype = 'float32')
    y_pred = preds * where #K.gather(preds,where)
    y_true =  y * where #K.gather(y,where)
    #print('s',y_pred.shape)
    loss = K.binary_crossentropy(y_true, y_pred)
    return loss

def semi_supervised_acc(y,preds):
    where = keras.backend.cast(y > -1 , dtype = 'float32')

    unlabled = K.sum(K.cast(y < 0, dtype='float32'))
    labled = K.sum(where)
    
    y_pred = preds * where #K.gather(preds,where)
    y_true =  y * where #K.gather(y,where)
    acc = keras.metrics.binary_accuracy(y_true, y_pred)
    
    return (acc * (labled + unlabled) - unlabled) / labled

    #acc = keras.backend.cast(y_true == y_pred , dtype = 'int32')
    #acc = K.mean(acc)
     