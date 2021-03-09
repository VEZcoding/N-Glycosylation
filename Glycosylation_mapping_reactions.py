# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:48:48 2021

@author: VEZcoding
"""
from pandas import DataFrame
enzymes = ['Man'] 
init_str = ('F%sM%sA%sG%sS%s')# %(0,i,j,z,n)) #začetni input
reactions_in = [] #list_glikanov
reactions_out = []
d = {}  #dictionary za dataframe
d['Glycan_in'] = []
d['---->'] = []
d['Glycan_out'] = []
mannose = [9,8,7,6,5,4]
c = 9 
for i in range(0,len(mannose)) :
    add_in =  init_str %(0,c,0,0,0)
    add_out =  init_str %(0,c-1,0,0,0)
    d['Glycan_in'].append(add_in)
    d['---->'].append('---->')
    d['Glycan_out'].append(add_out)
    reactions_in.append(add_in)
    reactions_out.append(add_out)
    c -= 1
    k = 0
    if c == 5 or c == 4:
        while k < 1: 
            add_in =  init_str %(0,c,k,0,0)
            add_out = init_str %(0,c,k+1,0,0)
            reactions_in.append(add_in)
            reactions_out.append(add_out)
            d['Glycan_in'].append(add_in)
            d['---->'].append('---->')
            d['Glycan_out'].append(add_out)
            k+=1
        z = 0    
        if k == 1:
            while z < 1: 
                add_in =  init_str %(0,c,k,z,0)
                add_out = init_str %(0,c,k,z+1,0)
                reactions_in.append(add_in)
                reactions_out.append(add_out)
                d['Glycan_in'].append(add_in)
                d['---->'].append('---->')
                d['Glycan_out'].append(add_out)
                z+=1 
            s = 0    
        if z == 1:
            while s < 1: 
                add_in =  init_str %(0,c,k,z,s)
                add_out = init_str %(0,c,k,z,s+1)
                reactions_in.append(add_in)
                reactions_out.append(add_out)
                d['Glycan_in'].append(add_in)
                d['---->'].append('---->')
                d['Glycan_out'].append(add_out)
                s+=1
    k = 0
    if c == 3:
        while k < 4:
            add_in =  init_str %(0,c,k,0,0)
            add_out = init_str %(0,c,k+1,0,0)
            reactions_in.append(add_in)
            reactions_out.append(add_out)
            d['Glycan_in'].append(add_in)
            d['---->'].append('---->')
            d['Glycan_out'].append(add_out)
            k+=1
            z = 0    
            if k != 0:
                while z < k: 
                    add_in =  init_str %(0,c,k,z,0)
                    add_out = init_str %(0,c,k,z+1,0)
                    reactions_in.append(add_in)
                    reactions_out.append(add_out)
                    d['Glycan_in'].append(add_in)
                    d['---->'].append('---->')
                    d['Glycan_out'].append(add_out)
                    z+=1 
                    s = 0    
                    if z != 0:
                        while s < z:
                            add_in =  init_str %(0,c,k,z,s)
                            add_out = init_str %(0,c,k,z,s+1)
                            reactions_in.append(add_in)
                            reactions_out.append(add_out)
                            d['Glycan_in'].append(add_in)
                            d['---->'].append('---->')
                            d['Glycan_out'].append(add_out)
                            s+=1 
            
        
df = DataFrame.from_dict(d)   #ustvari dataframe  
print(df) #izris rezultata 