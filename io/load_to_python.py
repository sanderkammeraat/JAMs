#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 13:35:38 2026

@author: kammeraat
"""

import h5py as h5

import os
import numpy as np

import matplotlib.pyplot as plt
#%%



save_folder_path = "/Users/kammeraat/mounting/pi-henkes/kammeraat/pairAN/sim_for_hdf5reader/simdata/"



raw_data = h5.File(os.path.join(save_folder_path, "raw_data.h5"), 'r')


#%%

frames = raw_data["frames"]



integration_info = raw_data["integration_info"]

integration_info.keys()
integration_info["Tsave"]


x1 = np.array(frames["1"]["x"])

integration_info["dt"]


#%%



raw_data["system"].keys()

raw_data["sizes"]







#%%


print(raw_data["frames"]["1"].keys())

#%% Suppose we want to collect the particle coordinates over time:
    
    
Nframes = len(raw_data['frames'].keys())
Np = len(raw_data['frames']['1']['x'])
    
x = np.zeros( shape=(Nframes, Np))

y = np.zeros( shape=(Nframes, Np))

ids = np.zeros( shape=(Nframes, Np))
types = np.zeros( shape=(Nframes, Np))

frame_numbers = np.arange(1,Nframes+1)



#%%




for i, frame_int in enumerate(frame_numbers):
    
    #As we can observe by printing the keys of raw_data['frames'] in order,
    #the key does not match up with its order (i.e. i!=int(frame))
    frame = str(frame_int)
    print(frame)
    
    #We correct for this by indexing the pre-allocated arrays with the actual frame key.
    #Note that frame 1 should be at index 0.
    
    x[i,:] = raw_data['frames'][frame]['x']
    y[i,:] = raw_data['frames'][frame]['y']
    
    ids[i,:] = raw_data['frames'][frame]['id']
    
    types[i,:] = raw_data['frames'][frame]['type']
    

    
#%% Let's plot the results to see if it makes sense
i = 200

plt.figure()

plt.scatter(x[i,:],y[i,:],c=ids[i,:])

plt.gca().set_aspect("equal")

plt.show()   
#%%


raw_data["system"].keys()


system =raw_data["system"]


force_params = {}

for forcetype in system["forces"]:
    
    
    if forcetype not in force_params:
        force_params[forcetype] = {}
    
    for force in  system["forces"][forcetype]:
        
        if force not in force_params[forcetype]:
            force_params[forcetype][force] = {}
        
        
        
        
        for field in system["forces"][forcetype][force]:
            

            force_params[forcetype][force][field] = np.array(system["forces"][forcetype][force][field])




#%%
raw_data.close()