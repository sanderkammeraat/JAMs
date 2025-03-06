#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#Here we show how the raw_data.jld2 hdf5 files can be loaed in python. 
#Besides the standard packages,
import numpy as np
import matplotlib.pyplot as plt


#we will need the following python package to load in hdf5 files, which in our case just happens to have a .jld2  file extension.
import h5py

#%%

#Let's load in the file
raw_data = h5py.File("/Users/kammeraat/JAMs_tutorial/raw_data.jld2")


#And show what it contains
print(raw_data.keys())

#As in the Julia tutorial, let's load in the x-coordinates of  frame number 20 in a numpy array.
x_20 = np.array(raw_data["frames"]["20"]["x"])

print(x_20)


#Alternatively, we can extract it directly by slicing with empy brackets
x_20_slice = raw_data["frames"]["20"]["x"][()]
print(x_20_slice)



#And likewise the time of the this frame
t_20 = raw_data["frames"]["20"]["t"][()]
print(t_20)
print('\n')


#%% Suppose we want to collect the particle coordinates over time:
x = np.zeros( shape=(len(raw_data['frames'].keys()), len(raw_data['frames']['1']['x'])))

y = np.zeros( shape=(len(raw_data['frames'].keys()), len(raw_data['frames']['1']['y'])))

ids = np.zeros( shape=(len(raw_data['frames'].keys()), len(raw_data['frames']['1']['id'])))
types = np.zeros( shape=(len(raw_data['frames'].keys()), len(raw_data['frames']['1']['type'])))




for i, frame in enumerate(raw_data['frames'].keys()):
    
    #As we can observe by printing the keys of raw_data['frames'] in order,
    #the key does not match up with its order (i.e. i!=int(frame))
    print(frame)
    
    #We correct for this by indexing the pre-allocated arrays with the actual frame key.
    #Note that frame 1 should be at index 0.
    
    x[int(frame)-1,:] = raw_data['frames'][frame]['x']
    y[int(frame)-1,:] = raw_data['frames'][frame]['y']
    
    ids[int(frame)-1,:] = raw_data['frames'][frame]['id']
    
    types[int(frame)-1,:] = raw_data['frames'][frame]['type']
    
    
    

    
#%% Let's plot the results to see if it makes sense
i = 10

plt.figure()

plt.scatter(x[i,:],y[i,:],c=ids[i,:])

plt.show()   
    

