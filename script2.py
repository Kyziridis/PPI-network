#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:00:58 2018

@author: dead

DBDM - Assignment.4 Protein-Protein-Interaction
"""
import matplotlib.pyplot as plt
import plotly.plotly as py
import networkx as nx
import numpy as np
import csv
from tqdm import tqdm
import itertools as it
import os

os.chdir("/home/dead/Documents/DBDM/Assignment-4/DM_Assignment_02/ResourceFiles")
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
print('\nCalculating the weight for each protein-connection as the portion of the intersection of their functions  >_\n')
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

THRS = input("\nProvide a number for threshold:  ")

# Find the couples above THRESHOLD :  to be cancer candidates
candidate_inx = np.where(neighbors[:,1] > int(THRS))
non_cand_inx = np.where(neighbors[:,1] <= int(THRS))
cancer_candidates = neighbors[candidate_inx]
non_cancer_candidates = neighbors[non_cand_inx]



# first prediction testdn are predicted as cancers
# testdn = prediction of cancers according to threshold
testdn = []
for i in range(len(cancer_candidates)): 
    if cancer_candidates[i][0][1] not in cancer:
        testdn.append([cancer_candidates[i][0][1] , "cancer"] )    
# testdn_non = prediction of NON_cancer according to the threshold        
testdn_non = []
for i in range(len(non_cancer_candidates)): 
    if non_cancer_candidates[i][0][1] not in cancer:
        testdn_non.append([non_cancer_candidates[i][0][1] , "nonCancer"] )

# Just the names of the predicted cancer_candidates
names_of_dirneighbors = np.asarray(testdn)[:,0]

# Evaluation of first naiv prediction according to the 'common' proteins in test1.txt        
name_cancer = np.intersect1d(testdn, truth[:,0]) 
name_cancer = np.reshape(name_cancer , (len(name_cancer),1))
indx = np.where(name_cancer == truth[:,0])[1] # Common proteins from prediction with truth

# Evaluation of NON_cancer predictions
name_nonCancer = np.intersect1d(testdn_non , truth[:,0])
name_nonCancer = np.reshape(name_nonCancer , (len(name_nonCancer),1 ) )
indx_non = np.where(name_nonCancer == truth[:,0])[1] # Common proteins from prediction of NONcancers with truth

# Conscructing the evaluation elements: TruePositive - FalsePositive - FalseNegative
correct = np.where(truth[indx,1] == 'cancer') # TruePositive
fp = np.where(truth[indx,1] == 'nonCancer') # FalsePositive (pred cancer -> truth nonCancer)
fn = np.where(truth[indx_non,1] == 'cancer')# FalseNegative (pred NONcancer -> truth Cancer)

# Calculate the evaluation formulas: Precision - Recall and F-Measure
Precision = correct[0].shape[0]/(len(fp[0]) +correct[0].shape[0])  
Recall = len(correct[0])/(len(correct[0]) + len(fn[0]))
F_measure = 2*Precision*Recall/(Precision + Recall)


print("\nFirst Naiv Approach: Direct neighbors of cancer proteins                                 >_")
print('\nPrecision: ' + str(Precision) )
print('\nRecall: ' + str(Recall) + ' F-Measure: ' + str(F_measure) )
################################################################################################





#################################################################
########################################################
# Shortest Path Algorithm #
print("\n")
print("Second Approach: Finding shortest path from all cancer variables")
print('\nConstructing the shortest_path_list for all cancer proteins    >_\n')
#comb = list(it.combinations(range(len(cancer)),2))
#comb = list(it.combinations(range(names_of_dirneighbors.shape[0]),2))
srtlist = []
for i in tqdm(cancer):
    p = nx.shortest_path(G, source=i , weight=None)
    srtlist.append(p)


# Flatten the srtlist with sortest paths
flat_list = [item for sublist in srtlist for item in sublist]       

# Cleaning srtlist from already predicted cancer neihbors, leaving inside only new predictions
# Flatten the object dirneigh which contains cancers & direct_neihbors (first prediction)
flt = dirneigh.flatten()        
un = np.unique(flt) #unique proteins which are already predicted as cancers or already known as cancers


# We need the proteins in flat_list which are not in dirneigh (already predicted).
# This is our new prediction for cancer_candidates
pred_new = np.setdiff1d(np.unique(flat_list),un) # array of values in flat_list that are not in dirneigh.
pred_new = np.reshape(pred_new , (len(pred_new),1))

# Prediction of NON_cancer
pred_nC = np.setdiff1d(G.nodes , pred_new )
pred_nC = np.reshape(pred_nC , (len(pred_nC),1) )

# Make a vector repeating 'cancer' for predictions and 'nonCancer'
c = np.repeat('cancer' , len(pred_new))
c = np.reshape(c , (len(c) ,1 ))
non_c = np.repeat('nonCancer' , len(pred_nC))
non_c = np.reshape(non_c , (len(non_c),1) )

# Construct the prediction list with the flag 'cancer' and 'nonCancer' next to each protein
testdn_new = np.hstack((pred_new ,c ) )
testdn_new_non = np.hstack( (pred_nC , non_c) )

# Evaluation of the second prediction for Cancers
name1 = np.intersect1d(testdn_new[:,0], truth[:,0])
name1 = np.reshape(name1 , (len(name1),1))
indx1 = np.where(name1 == truth[:,0])[1]

# Evaluation of the second prediction for nonCancers
name1_non = np.intersect1d(testdn_new_non[:,0], truth[:,0])
name1_non = np.reshape(name1_non , (len(name1_non),1))
indx1_non = np.where(name1_non == truth[:,0])[1]

# Conscructing the evaluation elements: TruePositive - FalsePositive - FalseNegative
correct1 = np.where(truth[indx1,1] == 'cancer') # TruePositive
fp1 = np.where(truth[indx1,1] == 'nonCancer') # FalsePositive (pred cancer -> truth nonCancer)
fn1 = np.where(truth[indx1_non,1] == 'cancer')# FalseNegative (pred NONcancer -> truth Cancer)

# Calculate the evaluation formulas: Precision - Recall and F-Measure
Precision = correct1[0].shape[0]/truth[indx1].shape[0] 
Recall = len(correct1[0])/(len(correct1[0]) + len(fn1[0]))
F_measure = 2*Precision*Recall/(Precision + Recall)



print("\nCalculating the weights for each connection through shortest path and the sum of weights")
# Make the list of all weights from shortest paths
telikol = []
for can in tqdm(range(len(srtlist))):
    
    adreas = []
    for i,j in list(srtlist[can].items()):
        if len(j)>1:
            metritis = 0 
            for k in range(len(j)-1):
                
                
                a = list(G[j[k]][j[k+1]].values())[0]
                
                metritis = metritis + a
            adreas.append([i , metritis])
    #teliko = np.append(teliko , adreas)
    telikol.append(adreas)
    
    
    
    
# Sort all lists inside telikol_list
#add_weigths = []    
r1 = np.zeros((len(telikol[0])))
#r1 = np.reshape(r1 , (len(r1) ,1))
for i in tqdm(range(len(telikol))):
    telikol[i].sort()
    r = list(map(int,np.asarray(telikol[i])[:,1])) # weights as integers from the list
    r1 = r1 + np.asarray(r)
    
# just the name of the proteins sorted    
prots  = np.reshape(np.asarray(telikol[0])[:,0] , (len(r1),1))
r1 = np.reshape(r1 , (len(r1),1))  
 
add_weigths = np.hstack((prots,r1 ) )
    


a = list(map(float,add_weigths[:,1]))
flag = input("\nPlot the histogram of sum of weights? [y/n]:  ")
if flag == 'y':
    plt.figure()
    plt.hist(a)
    plt.show()    
THRS2 = input('\nProvide a number for threshold e.g. 10000:  ')   

lista = []     
for i in range(len(a)):
    if a[i] > int(THRS2):
        lista.append(['cancer'])
    else:
        lista.append(['nonCancer'])
    
    
soupa = np.hstack((prots , lista))    

cancer_indx = np.where(soupa[:,1] == 'cancer')
non_cancer_indx = np.where(soupa[:,1] == 'nonCancer')      

names1_cancer = np.intersect1d(soupa[cancer_indx][:,0] , truth[:,0])  
names1_cancer = np.reshape(names1_cancer , (len(names1_cancer),1))  
indx_cancer_e = np.where(names1_cancer == truth[:,0])[1]
    
names1_non_cancer = np.intersect1d(soupa[non_cancer_indx][:,0] , truth[:,0])  
names1_non_cancer = np.reshape(names1_non_cancer , (len(names1_non_cancer),1))  
indx_non_cancer_e = np.where(names1_non_cancer == truth[:,0])[1]

    
tp = len(np.where(truth[indx_cancer_e,1] == 'cancer')[0]) 
fp = len(np.where(truth[indx_cancer_e,1] == 'nonCancer')[0]) 
fn = len(np.where(truth[indx_non_cancer_e,1] == 'cancer')[0])


Precision = tp/(fp + tp)
Recall = tp/(fn+tp)
F_measure = 2*Precision*Recall/(Precision + Recall)

print('\nPrecision: ' + str(Precision) )
print('\nRecall: ' + str(Recall) + ' F-Measure: ' + str(F_measure) )






















    
    
    
    