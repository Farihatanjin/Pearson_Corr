#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 00:33:47 2020

@author: farihatanjin
"""




import pandas as pd
import matplotlib.pyplot as plt
import csv
import seaborn as sns
import numpy as np



rows=dict()


#read from csv file

data= open('Sample.csv')
csv_data=csv.reader(data)

#select relevant genes from larger dataset
for row1 in csv_data:
    #for row2 in csv.reader(open('/Users/farihatanjin/Desktop/lab/Book4.csv')):
    for row2 in csv.reader(open('genelist.csv')):
        #check if gene in csv file is present in genelist
        if row1[0]==row2[0]:
            rows[row1[0]]=row1[1:]
            #append to dictionary if gene is found
            
#matrix from dictionary 
matrix=pd.DataFrame(data=rows)

matrix = matrix.astype(float)

#correlation of matrix
df=matrix.corr(method='pearson')
print(df)

  


df.to_csv("correlogram.csv")


#only show lower triangular heatmap
mask = np.zeros_like(df, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True


f, ax = plt.subplots(figsize=(25, 20)) 

heatmap = sns.heatmap(df, 
                      xticklabels=1,
                      yticklabels=1,
                      mask = mask,
                      square = True,
                      linewidths = .5,
                      cmap='PiYG',
                      cbar_kws={"shrink": .4, 'ticks' : [-1, -.5, 0, 0.5, 1]},
                      vmin = -0.5, 
                      vmax = 1,
                      annot = False,
                      )
                  
ax.set_yticklabels(df, rotation=0)
ax.set_xticklabels(df.columns, fontsize=9)
plt.xticks(fontsize=9, rotation=90)
plt.yticks(fontsize=9)
sns.set_style({'xtick.bottom': True}, {'ytick.left': True})

figure = heatmap.get_figure()    
figure.savefig('LB.png', dpi=400)


#positive correlations   
pos = np.where(df > 0.21)
pos = [(df.index[x], df.columns[y]) for x, y in zip(*pos)
                                        if x != y and x < y]
with open("positive_corr.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(pos)
    
#negative correlations  
neg = np.where(df < -0.05)
neg = [(df.index[x], df.columns[y]) for x, y in zip(*neg)
                                        if x != y and x < y]
with open("negative_corr.csv", "w", newline="") as g:
    writer = csv.writer(g)
    writer.writerows(neg) 

#only look at lower tringular data   
def get_redundant_pairs(df):
    duplicates = set()
   
    cols = df.columns
    for i in range(0, df.shape[1]):
        for j in range(0, i+1):
            duplicates.add((cols[i], cols[j]))
  
    return duplicates

def get_top_correlations(df, n):
    au_corr = df.unstack()
    dupl = get_redundant_pairs(df)
    au_corr = au_corr.drop(labels=dupl).sort_values(ascending=False).drop_duplicates()
    return au_corr[0:n]

def get_low_correlations(df, n):
    au_corr = df.unstack()
    dupl = get_redundant_pairs(df)
    au_corr = au_corr.drop(labels=dupl).sort_values(ascending=True).drop_duplicates()
    return au_corr[0:n]


get_top_correlations(df, 20).to_csv("top_corr.csv")

get_low_correlations(df, 20).to_csv("low_corr.csv")
  
