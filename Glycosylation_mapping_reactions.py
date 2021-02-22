# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:48:48 2021

@author: VEZcoding
"""
from pandas import DataFrame
enzymes = ['Man']   # DOPOLNI 
init_str = ('F%sM%sA%sG%sS%s')# %(0,i,j,z,n))
reactions = []
mannose = [9,8,7,6,5,4,3]
for i in mannose :
    add =  (init_str %(0,i,0,0,0))
    reactions.append(add)
    if i == 5 or i == 4: 
        for j in range(1,2):
            add = (init_str %(0,i,j,0,0))
            reactions.append(add)
            for z in range (1,2):
                if z > j: 
                    continue
                else: 
                    add = (init_str %(0,i,j,z,0))
                    reactions.append(add)
                for n in range (1,2):
                      if n > z: 
                          continue
                      else: 
                          add = (init_str %(0,i,j,z,n))
                          reactions.append(add)
    if i == 3:
        for j in range(1,5):
            add = (init_str %(0,i,j,0,0))
            reactions.append(add)
            for z in range (1,5):
                if z > j: 
                    continue
                else: 
                    add = (init_str %(0,i,j,z,0))
                    reactions.append(add)
                for n in range (1,5):
                      if n > z: 
                          continue
                      else: 
                          add = (init_str %(0,i,j,z,n))
                          reactions.append(add)
print(reactions)
d = {}
d['Glycan_in'] = []
d['---->'] = []
d['Glycan_out'] = []
c = 1
for i in range(0, len(reactions)+1):    
    d['Glycan_in'].append(reactions[i])
    d['---->'].append('--->')
    d['Glycan_out'].append(reactions[c])
    c += 1
    
    
df = DataFrame.from_dict(d)    
print(df)