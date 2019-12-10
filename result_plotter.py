
import pandas as pd
import matplotlib.pyplot as plt
import progressbar
import numpy as np

def get_vessel_points(vessel_data):
    name_coords = dict()
    for i in progressbar.progressbar(range(len(vessel_data))):
        point = np.array((vessel_data.iloc[i]['x'],vessel_data.iloc[i]['y'],vessel_data.iloc[i]['z']))
        name =  str(vessel_data.iloc[i]['nodeid']).strip()
        
        if name in name_coords.keys():
            #print(name)
            name_coords[name].append(point)
        else:
            name_coords[name] = [point]
            
    for key in progressbar.progressbar(name_coords.keys()):
        name_coords[key] = np.array(name_coords[key])
        
    return name_coords


def get_tumor_points(tumor_data):
    tumor = []
    for i in range(len(tumor_data)):
        point = np.array((tumor_data.iloc[i][0],tumor_data.iloc[i][1],tumor_data.iloc[i][2] +25))
        tumor.append(point)
    tumor = np.array(tumor)
    return tumor






lung_name = "lung01_041" 


results = pd.read_csv('predictions_'+lung_name+'.csv')
name_index = {str(results.id[i]).strip():i for i in range(len(results))}
name_pred = {str(results.id[i]).strip():results['0'][i] for i in range(len(results))}

vessel_data = pd.read_csv('D:\\work\Classes\\Tufts\\MLonGraphs\\Project\\data\\'+lung_name+'\\files\\'+lung_name+'_coords_scaled.csv', sep = ',')
tumor_data = pd.read_csv('D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\'+lung_name+'\\files\\'+lung_name+'_tumor_coords.csv', sep = ',')

features = pd.read_csv('D:\work\Classes\Tufts\MLonGraphs\Project\data\lung_adjmatrices_featurematrices\\'+lung_name+'_feature_matrix.csv')


tumor = get_tumor_points(tumor_data)
name_coords = get_vessel_points(vessel_data)

#in_tumor = get_vessels_in_tumor(tumor,name_coords) #to be thrown away because false reconstructions inside solid mass
#near_tumor = get_vessels_near_tumor(tumor,name_coords) # to be labeled diseased (need to remove in_tumor first)
name_labels = {str(features.id[i]).strip():features.labels[i] for i in range(len(features))}

#print(name_coords)
colors = ['blue','red']

from mpl_toolkits.mplot3d import axes3d, Axes3D
fig = plt.figure(1)
ax = Axes3D(fig)

ax.plot3D(tumor.T[0], tumor.T[1], tumor.T[2], 'gray')
for key in name_coords.keys():
    if not len(name_coords[key]) == 0:
        #print(name_coords[key])
        if key in name_pred.keys():
            if name_pred[key] == 1:
                
                ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'red')#colors[name_label[key]])
            if name_pred[key] == 0:
                True
                ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'blue')
            if name_labels[key] == -1:
                True
                ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'green')#

fig = plt.figure(2)
ax = Axes3D(fig)

ax.plot3D(tumor.T[0], tumor.T[1], tumor.T[2], 'gray')
for key in name_coords.keys():
    if not len(name_coords[key]) == 0:
        #print(name_coords[key])
        if key in name_labels.keys():
            if name_labels[key] == 1:
                
                ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'red')#colors[name_label[key]])
            if name_labels[key] == 0 :
                True
                ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'blue')
            if name_labels[key] == -1:
                True
                ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'green')#colors[name_label[key]])

   
fig = plt.figure(3)
ax = Axes3D(fig)

ax.plot3D(tumor.T[0], tumor.T[1], tumor.T[2], 'gray')
for key in name_coords.keys():
    if not len(name_coords[key]) == 0:
        #print(name_coords[key])
        if key in name_pred.keys():
            if name_pred[key] == 1:
                
                ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'red')#colors[name_label[key]])
            if name_pred[key] == 0:
                True
                ax.plot3D(name_coords[key].T[0],name_coords[key].T[1],name_coords[key].T[2], color = 'blue')
           
           

plt.show()




