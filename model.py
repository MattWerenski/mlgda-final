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
from MyGraphConv import *
from myloss import *
from tqdm import tqdm 
class SimpleModel:
    def __init__(self,MAX_GRAPH_SIZE,F):
        self.MAX_GRAPH_SIZE = MAX_GRAPH_SIZE
        
        A = Input((MAX_GRAPH_SIZE, ),) #sparse=True)
        X = Input((F,))
        
        #h1 = GraphConv(100)([X,A])
        
        h2 = GraphSageConv(20)([X,A])
        h2 = GraphSageConv(20)([X,A])
        
        
       

        d = Dense(64, activation = 'relu')(h2)
        y = Dense(1, activation = 'sigmoid')(d)
    
        self.model = Model(inputs= [A,X], outputs = y)

        learning_rate = 0.01
        optimizer = keras.optimizers.Adam()
        self.model.compile(optimizer=optimizer,
              loss=semi_supervised,
              metrics=[semi_supervised_acc])
        return 
        
    def fit(self,A,X, labels):
        print(A.shape,X.shape,labels.shape)
        self.model.fit([A,X], labels,
                epochs=10,
                batch_size=self.MAX_GRAPH_SIZE,
                validation_split = 0.0, verbose = 1, shuffle = False)
        return
    
    def fit_many(self,data):
        #print(A.shape,X.shape,labels.shape)
        for i in tqdm(range(10)):
            for datum in data:
                A,X, labels = datum
                #print(A.shape)
                self.model.fit([A,X], labels,
                        epochs=5,
                        batch_size=self.MAX_GRAPH_SIZE,
                        validation_split = 0.0, verbose = 0, shuffle = False)
        return

    def predict(self, A, X):
        """
        Predict labels for a test set
        """
        
        return self.model.predict([A,X], batch_size= self.MAX_GRAPH_SIZE, verbose=0)

    def evaluate(self, A, X, labels):
        """
        Evaluate the model on a test set
        """
        
       # adj = nx.to_numpy_array(A)
       
        score = self.model.evaluate([A,X], labels, batch_size=self.MAX_GRAPH_SIZE)

        return score



