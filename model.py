import tensorflow as tf
from tensorflow import keras
import spektral
import keras
from spektral.layers import GraphSageConv, GraphConv, GraphAttention
from keras.layers import Input, Dense, Dropout, Concatenate, Add
from keras.models import Model
from keras.optimizers import Adam
from keras.layers import Layer
from keras.activations import relu
from sklearn import metrics
#from MyGraphConv import *
from myloss import *
from tqdm import tqdm 
from sklearn.metrics import confusion_matrix, f1_score
class SimpleModel:
    def __init__(self,MAX_GRAPH_SIZE,F):
        self.MAX_GRAPH_SIZE = MAX_GRAPH_SIZE
        
        A = Input((MAX_GRAPH_SIZE, ),) #sparse=True)
        X = Input((F,))
        
        #h1 = GraphConv(50, activation='relu')([X,A])
        
        #h1 = GraphConv(20, activation='relu')([h1,A])
        #h1 = GraphConv(20, activation='relu')([h1,A])
        #h1 = GraphConv(20, activation='relu')([h1,A])
        #h1 = GraphConv(20, activation='relu')([h1,A])
        #h1 = GraphConv(10)([h1,A])

        h2 = GraphSageConv(100, activation='relu')([X,A])


        #h2 = GraphSageConv(50, activation='relu')([h2,A])
        
        h2 = GraphSageConv(100, activation='relu')([h2,A])
        h2 = GraphSageConv(100, activation='relu')([h2,A])
        h2 = GraphSageConv(50, activation='relu')([h2,A])
        h2 = GraphSageConv(50, activation='relu')([h2,A])
        h2 = GraphSageConv(50, activation='relu')([h2,A])
        h2 = GraphSageConv(50, activation='relu')([h2,A])
        h2 = GraphSageConv(50, activation='relu')([h2,A])
        #h2 = GraphSageConv(10, activation='relu')([h2,A])
        
        #drop = Dropout(0.5)(X)
        #h2 = GraphSageConv(20, activation='relu')([drop,A])
        
        #h2 = GraphSageConv(20, activation='relu')([h2,A])
        
        #drop = Dropout(0.5)(h2)
        #h2 = GraphSageConv(50,activation='relu')([h2,A])
        

        #drop = Dropout(0.5)(h2)
        #h2 = GraphSageConv(50)([h2,A])
        #h3 = GraphAttention(20)([h2,A])
        #drop = Dropout(0.5)(h3)
        #h3 = GraphAttention(20)([drop,A])
        #drop = Dropout(0.5)(h3)
        #h3 = GraphAttention(20)([drop,A])
        

        #d = Dense(10,activation='relu')(h2)
        #conc = Concatenate()([h2,d])

        y = Dense(1, activation = 'sigmoid')(h2)
    
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
                epochs=50,
                batch_size=self.MAX_GRAPH_SIZE,
                validation_split = 0.0, verbose = 1, shuffle = False, class_weight={-1:0,0:1,1:4})
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
        preds=  self.predict(A,X)

        semi_preds = preds[np.where(labels>-1)]
        int_preds = np.round(preds)
        semi_labels = labels[np.where(labels>-1)]
        semi_intpreds = int_preds[np.where(labels>-1)]
        
        
        conf = confusion_matrix(semi_labels,semi_intpreds).ravel()
        fpr, tpr, thresholds = metrics.roc_curve(semi_labels, semi_preds.reshape(-1), pos_label=1)
        auc = metrics.auc(fpr, tpr)
        f1 = f1_score(semi_labels, semi_intpreds, zero_division=1)
        return score, conf, auc,f1

class TwoAdjModel:
    def __init__(self,MAX_GRAPH_SIZE,F):
        self.MAX_GRAPH_SIZE = MAX_GRAPH_SIZE
        
        Av = Input((MAX_GRAPH_SIZE, ),) #sparse=True)
        As = Input((MAX_GRAPH_SIZE, ),)
        X = Input((F,))
        
        # GraphSage layers for the vascular network
        '''
        hv = GraphSageConv(100, activation='relu')([X,Av])
        hv = GraphSageConv(100, activation='relu')([hv,Av])
        hv = Dropout(0.5)(hv)
        hv = GraphSageConv(50, activation='relu')([hv,Av])
        hv = GraphSageConv(50, activation='relu')([hv,Av])
        hv = Dropout(0.5)(hv)
        hv = GraphSageConv(50, activation='relu')([hv,Av])
        hv = Dropout(0.5)(hv)

        # GraphSage layers for the spatial network
        hs = GraphSageConv(100, activation='relu')([X,As])
        hs = GraphSageConv(100, activation='relu')([hs,As])
        hs = Dropout(0.5)(hs)
        hs = GraphSageConv(50, activation='relu')([hs,As])
        hs = GraphSageConv(50, activation='relu')([hs,As])
        hs = Dropout(0.5)(hs)
        hs = GraphSageConv(50, activation='relu')([hs,As])
        hs = Dropout(0.5)(hs)
        

        con = Concatenate()([hv, hs])
        '''

        hv = GraphSageConv(100, activation='relu')([X,Av])
        hs = GraphSageConv(100, activation='relu')([X,As])
        s = Add()([hv,hs])
        hv = GraphSageConv(100, activation='relu')([s,Av])
        hs = GraphSageConv(100, activation='relu')([s,As])
        s = Add()([hv,hs])
        hv = GraphSageConv(50, activation='relu')([s,Av])
        hs = GraphSageConv(50, activation='relu')([s,As])
        s = Add()([hv,hs])
        hv = GraphSageConv(50, activation='relu')([s,Av])
        hs = GraphSageConv(50, activation='relu')([s,As])
        s = Add()([hv,hs])
        hv = GraphSageConv(50, activation='relu')([s,Av])
        hs = GraphSageConv(50, activation='relu')([s,As])
        s = Add()([hv,hs])

        y = Dense(1, activation = 'sigmoid')(s)
    
        self.model = Model(inputs= [Av, As, X], outputs = y)

        learning_rate = 0.01
        optimizer = keras.optimizers.Adam()
        self.model.compile(optimizer=optimizer,
              loss=semi_supervised,
              metrics=[semi_supervised_acc])
        return 
        
    def fit(self,Av, As, X, labels):
        print(Av.shape, As.shape, X.shape,labels.shape)
        self.model.fit([Av, As, X], labels,
                epochs=40,
                batch_size=self.MAX_GRAPH_SIZE,
                validation_split = 0.0, verbose = 1, shuffle = False, class_weight={-1:0,0:1,1:4})
        return
    

    def predict(self, Av, As, X):
        """
        Predict labels for a test set
        """
        
        return self.model.predict([Av, As, X], batch_size= self.MAX_GRAPH_SIZE, verbose=0)

    def evaluate(self, Av, As, X, labels):
        """
        Evaluate the model on a test set
        """
        
       # adj = nx.to_numpy_array(A)
       
        score = self.model.evaluate([Av, As, X], labels, batch_size=self.MAX_GRAPH_SIZE)
        preds=  self.predict(Av, As, X)

        semi_preds = preds[np.where(labels>-1)]
        int_preds = np.round(preds)
        semi_labels = labels[np.where(labels>-1)]
        semi_intpreds = int_preds[np.where(labels>-1)]
        
        
        conf = confusion_matrix(semi_labels,semi_intpreds).ravel()
        fpr, tpr, thresholds = metrics.roc_curve(semi_labels, semi_preds.reshape(-1), pos_label=1)
        auc = metrics.auc(fpr, tpr)
        f1 = f1_score(semi_labels, semi_intpreds, zero_division=1)
        return score, conf, auc,f1


class DenseModel:
    def __init__(self,MAX_GRAPH_SIZE,F):
        self.MAX_GRAPH_SIZE = MAX_GRAPH_SIZE
        
        A = Input((MAX_GRAPH_SIZE, ),) #sparse=True)
        X = Input((F,))
        
        
        

        d = Dense(60, activation = 'relu')(X)
        d = Dense(20, activation = 'relu')(d)
        d = Dense(20, activation = 'relu')(d)
        d = Dense(20, activation = 'relu')(d)
        d = Dense(20, activation = 'relu')(d)

        
        d = Dense(20, activation = 'relu')(d)
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
                epochs=100,
                batch_size=self.MAX_GRAPH_SIZE,
                validation_split = 0.0, verbose = 1, shuffle = False, class_weight={-1:0,0:1,1:3})
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
        preds=  self.predict(A,X)

        int_preds = np.round(preds)
        semi_labels = labels[np.where(labels>-1)]
        semi_intpreds = int_preds[np.where(labels>-1)]
    

        conf = confusion_matrix(semi_labels,semi_intpreds).ravel()

        return score, conf




