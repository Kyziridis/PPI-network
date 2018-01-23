#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:00:58 2018

@author: dead

DBDM - Assignment.4 Protein-Protein-Interaction
"""
import pickle
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import csv
from tqdm import tqdm
from collections import defaultdict
from collections import Counter
import plotly.plotly as py
import random
import itertools as it

# open the file humanPPI.txt and make it undirected
G1=nx.read_edgelist('humanPPI.txt', delimiter=',',create_using=nx.DiGraph())
Gu=G1.to_undirected()
f =Gu.edges(Gu)
f =np.array([(x) for x in f]) # f object contains the protein undirected connections


# Import files
cancer = np.loadtxt("Cancer.txt",dtype=str)
function = np.genfromtxt("Functions.txt" , dtype=str , delimiter=',' , usecols=(0,1))
truth = np.genfromtxt("Test1.txt" , dtype=str , delimiter=',' , usecols=(0,1))
pred = np.genfromtxt("Test2.txt" , dtype=str , delimiter=',' , usecols=(0,1))
#############################################################################

# Calculate the weights for each P-P connection according to the intersection of their functions
print('\nCalculating the weight for each protein-connection\n')
l = []
for i in tqdm(range(len(f))) :
    l.append(len(np.intersect1d(function[f[i][0] == function[:,0]][:,1] , function[f[i][1] == function[:,0]][:,1])) )
    

# Transpose list l to an array and construct the newdb as an array matrix with the hole network
ll = np.array(l)
ll = np.reshape(ll , (f.shape[0],1))
newdb = np.hstack((f,ll))

# Then extract the newdb as a .csv file in order to work with.    
with open('dbdm.csv', 'w',encoding="utf-8") as csv_file:
    writer = csv.writer(csv_file, delimiter='\t',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for j in newdb:
        writer.writerow(j)


# Read the new data base 
G=nx.read_edgelist('dbdm.csv', delimiter='\t',create_using=nx.Graph(),data=[('weight',int)])
nx.info(G)

# Find dierct neihbors only for the cancer-proteins
dirneigh = []
for i in cancer:
    Out = G.neighbors(i)
    for j in Out:
        dirneigh.append([i,j])

dirneigh = np.asarray(dirneigh)
# MAking cancer nodes and weights
neighbors = np.empty(0)
for i in np.asarray(dirneigh):
    w = G.get_edge_data(*tuple(i))
    neighbors=np.append(neighbors,[i,list(w.values())[0]] )
    
# Object neighbors is an array of cancer proteins and their direct neihbors with their weight.
neighbors = np.reshape(neighbors , (len(dirneigh) , 2))

# Average weight = 3.41
np.mean(neighbors[:,1])

# Find the couples above THRESHOLD : 
candidate_inx = np.where(neighbors[:,1] > 8)
cancer_candidates = neighbors[candidate_inx]

# testdn = prediction of cancers according to threshold
testdn = []
for i in range(len(cancer_candidates)): 
    if cancer_candidates[i][0][1] not in cancer:
        testdn.append([cancer_candidates[i][0][1] , "cancer"] )    
# first prediction testdn are predict_cancers        


# Evaluation of first naiv prediction according to the 'common' proteins in test1.txt
name = np.intersect1d(testdn, truth[:,0])
name = np.reshape(name , (len(name),1))
indx = np.where(name == truth[:,0])[1]
truth[indx]

correct = np.where(truth[indx,1] == 'cancer')
print('\nAccuracy in '+str(truth[indx].shape[0])+' sample: ' + str(correct[0].shape[0]/truth[indx].shape[0]))
################################################################################################



########################################################
# Shortest Path Algorithm #

print('\nConstructing the shortest path list above 4 neighbors    >_\n')
comb = list(it.combinations(range(len(cancer)),2))
srtlist = []
for i,j in tqdm(comb):
    p = nx.shortest_path(G, source=cancer[i],target=cancer[j] , weight=None)
    if len(p) > 4:
        srtlist.append(p)

sort = np.array(srtlist)

# Flatten the srtlist 
flat_list = [item for sublist in srtlist for item in sublist]       

# Cleaning srtlist from already predicted cancer neihbors, leaving inside only new predictions
# Flatten the object dirneigh which contains cancers-direct_neihbors (first prediction)
flt = dirneigh.flatten()        
un = np.unique(flt) #unique proteins which are already predicted as cancers or already known as cancers

# Sorted 1D array of values in ar1 that are not in ar2.
# We need the proteins in flat_list which are not in dirneigh.
# This is our new prediction about for cancer_candidates
pred_new = np.setdiff1d(np.unique(flat_list),un)
pred_new = np.reshape(pred_new , (len(pred_new),1))

# Make a vector repeating 'cancer' for predictions
c = np.repeat('cancer' , len(pred_new))
c = np.reshape(c , (len(c) ,1 ))

# Construct the prediction list with the flag 'cancer' next to each protein
testdn_new = np.hstack((pred_new ,c ) )

# Evaluation of the second prediction
name1 = np.intersect1d(testdn_new[:,0], truth[:,0])
name1 = np.reshape(name1 , (len(name1),1))
indx1 = np.where(name1 == truth[:,0])[1]
truth[indx1]

correct1 = np.where(truth[indx1,1] == 'cancer')
print('Accuracy in '+str(truth[indx1].shape[0])+' sample: ' + str(correct[0].shape[0]/truth[indx1].shape[0]))










    







# Run sortest path
# =============================================================================
# sortest1 = {}
# for i in tqdm(cancer):    
#     p=nx.shortest_path(G,source=i)
#     sortest1[i].append(p)
# =============================================================================
    
    
    
 
    
    
    
    
    
    
    
    
    
    
    
    