import tensorflow as tf
from tensorflow import keras 
import pandas as pd
import numpy as np
from model import SimpleModel
from scipy.sparse import lil_matrix
from load_data import *

A,X,y = load_data()

model = SimpleModel(A.shape[0],X.shape[-1])

model.fit(A,X,y)



score = model.evaluate(A,X,y)

print(score)

