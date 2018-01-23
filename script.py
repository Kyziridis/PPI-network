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
        dirneigh.append((i,j))


# MAking cancer nodes and weights
x = tuple(dirneigh)
neighbors = np.empty(0)
edges = []
for i in x:
    w = G.get_edge_data(*i)
    neighbors=np.append(neighbors,[i,list(w.values())[0] ] )
    edges.append(i + (list(w.values())[0],))

neighbors = np.reshape(neighbors , (len(x) , 2))

# Object neighbors is an array of cancer proteins and their direct neihbors with their weight.





# Run sortest path
sortest1 = {}
for i in tqdm(cancer):    
    p=nx.shortest_path(G,source=i)
    sortest1[i].append(p)
    
    
    
 
    
    
    
    
    
    
    
    
    
    
    
    