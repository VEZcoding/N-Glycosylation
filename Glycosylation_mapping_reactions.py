# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:48:48 2021

@author: VEZcoding
"""
from pandas import DataFrame
import re
enzymes = ['Man'] 
init_str = ('F%sM%sA%sG%sS%s')# %(0,i,j,z,n)) #začetni input
reactions = [] #list_glikanov
d = {}  #dictionary za dataframe
d['Glycan_in'] = []
d['---->'] = []
d['Glycan_out'] = []
mannose = [9,8,7,6,5,4,3] 
for i in mannose :
    add =  init_str %(0,i,0,0,0)   
    reactions.append(add)
    #odvzemaj manoze in dopolni z eno vejo legokock
    if i == 5 or i == 4: 
        for j in range(0,2):
            #dodajaj A blocke
            add = init_str %(0,i,j,0,0)           
            reactions.append(add)
            #dodajaj G blocke
            for z in range (0,2):
                if z > j  and z: 
                    continue
                else: 
                    add = init_str %(0,i,j,z,0)                   
                    reactions.append(add)
                #dodajaj S blocke
                for n in range (0,2):
                      if n > z: 
                          continue
                      else: 
                          add = init_str %(0,i,j,z,n) 
                          reactions.append(add)
    if i == 3:
        #dodaj A
        for j in range(1,5):
            add = init_str %(0,i,j,0,0)
            reactions.append(add)
            #dodaj G
            for z in range (1,5):
                if z > j: 
                    continue
                else: 
                    add = init_str %(0,i,j,z,0)
                    reactions.append(add)
                #dodaj S
                for n in range (1,5):
                      if n > z: 
                          continue
                      else: 
                          add = init_str %(0,i,j,z,n)
                          reactions.append(add)
                          

reactions_1 = reactions[:-1]  #zadnj produkt != reaktant
reactions_2 = reactions[1:]   #začetni reaktant != produkt

for i in range(0, len(reactions)-1):
    r_1 = re.findall(r'\d+',reactions[i])
    r_2 = re.findall(r'\d+',reactions_2[i])
    print (r_1)
    if r_1[4] > r_2[4]:
        continue
    if r_1 == r_2: 
        continue
    else: 
        d['Glycan_in'].append(reactions_1[i])
        d['---->'].append('--->')
        d['Glycan_out'].append(reactions_2[i]) 
df = DataFrame.from_dict(d)   #ustvari dataframe  
print(df) #izris rezultata  