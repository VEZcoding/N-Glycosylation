# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:48:48 2021

@author: VEZcoding
"""
"""
Starting input of Unglycosylated product will be F0M0A0G0S0. Where each letter presents a glycosylation(building) block and number that follows presents the amount of
the blocks on the glycan - F=Fucose, M=Manose, A=GlucNAc, G = Galactose, S = Sialic Acid. For example F0M3A2G1S1 glycan is afucosylated(has no fucose), 3 mannose
2 GlucNACs and 1 Galactose and 1 Sialic Acid.  
List of enzymes present: 
    - ManI
    - GalT
    - GnT
    - SiaT
For more accurate N-Glycoslyation pathway we would have acknowledge more enzmyes. Like GnTIII, or GnTIV/GnTV, which can be incorporated, 
but in my opinion it is not necessary for a basic understanding of this pathway.     
    """

from pandas import DataFrame
enzymes = ['ManI','GalT','GnT','SiaT'] 
init_str = ('F%sM%sA%sG%sS%s') #starting input which will be F0M0A0G0S0
reactions_in = [] #list of gylcans
reactions_out = [] #list of glycan
d = {}  #dictionary for creating dataframe

#different keys in a dataframe
d['Gene'] = []
d['Enzyme'] = []
d['Glycan_in'] = []
d['Donor']  = []
d['---->'] = []
d['Glycan_out'] = []
d['Residue']  = []

mannose = [9,8,7,6,5,4] #number of mannose we would like to have
m = mannose[0]
for i in range(0,len(mannose)) : #extract Mannose
    add_in =  init_str %(0,m,0,0,0)
    add_out =  init_str %(0,m-1,0,0,0)
    d['Gene'].append('ManI')
    d['Enzyme'].append('ManI')
    d['Glycan_in'].append(add_in)
    d['Donor'].append('H2O')
    d['---->'].append('---->')
    d['Glycan_out'].append(add_out)
    d['Residue'].append('Man')
    reactions_in.append(add_in)
    reactions_out.append(add_out)
    m -= 1
    a = 0
    if m == 5 or m == 4: #append GlucNAc
        while a < 1: 
            add_in =  init_str %(0,m,a,0,0)
            add_out = init_str %(0,m,a+1,0,0)
            reactions_in.append(add_in)
            reactions_out.append(add_out)
            d['Gene'].append('GnTI')
            d['Enzyme'].append('GnTI')
            d['Glycan_in'].append(add_in)
            d['Donor'].append('UDP_GlcNAc')
            d['---->'].append('---->')
            d['Glycan_out'].append(add_out)
            d['Residue'].append('UDP')           
            add_in =  init_str %(0,m,a+1,0,0)
            add_out = init_str %(0,m-1,a+1,0,0)
            d['Gene'].append('ManII')
            d['Enzyme'].append('ManII')
            d['Glycan_in'].append(add_in)
            d['Donor'].append('UDP_GlcNAc')
            d['---->'].append('---->')
            d['Glycan_out'].append(add_out)
            d['Residue'].append('UDP')
            a+=1
        g = 0    
        if a == 1:
            if g == 0: 
                add_in =  init_str %(1,m,a,g,0)
                add_out = init_str %(1,m,a,g,0)
                reactions_in.append(add_in)
                reactions_out.append(add_out)
                d['Gene'].append('FucT')
                d['Enzyme'].append('FucT')
                d['Glycan_in'].append(add_in)
                d['Donor'].append('GDP_Fuc')
                d['---->'].append('---->')
                d['Glycan_out'].append(add_out)
                d['Residue'].append('GDP')
            while g < 1: #append Galactose
                add_in =  init_str %(0,m,a,g,0)
                add_out = init_str %(0,m,a,g+1,0)
                reactions_in.append(add_in)
                reactions_out.append(add_out)
                d['Gene'].append('GalT')
                d['Enzyme'].append('GalT')
                d['Glycan_in'].append(add_in)
                d['Donor'].append('UDP_Gal')
                d['---->'].append('---->')
                d['Glycan_out'].append(add_out)
                d['Residue'].append('UDP')
                g+=1

        s = 0    
        if g == 1:
            while s < 1: #append Sialic acid
                add_in =  init_str %(0,m,a,g,s)
                add_out = init_str %(0,m,a,g,s+1)
                reactions_in.append(add_in)
                reactions_out.append(add_out)
                d['Gene'].append('SiaT')
                d['Enzyme'].append('SiaT')
                d['Glycan_in'].append(add_in)
                d['Donor'].append('CMP_Neu5Ac')
                d['---->'].append('---->')
                d['Glycan_out'].append(add_out)
                d['Residue'].append('CMP')
                s+=1
    a = 0
    if m == 3: #if all manose is removed
        while a < m+1: #append GluNac m+1 indicates that 2 GlcNAc append to 1 branch of mannose
            add_in =  init_str %(0,m,a,0,0)
            add_out = init_str %(0,m,a+1,0,0)
            reactions_in.append(add_in)
            reactions_out.append(add_out)
            d['Gene'].append('GnTI')
            d['Enzyme'].append('GnTI')
            d['Glycan_in'].append(add_in)
            d['Donor'].append('UDP_GlcNAc')
            d['---->'].append('---->')
            d['Glycan_out'].append(add_out)
            d['Residue'].append('UDP')
            a+=1
            g = 0    
            if a != 0:
                if g == 0: 
                    add_in =  init_str %(1,m,a,g,0)
                    add_out = init_str %(1,m,a,g,0)
                    reactions_in.append(add_in)
                    reactions_out.append(add_out)
                    d['Gene'].append('FucT')
                    d['Enzyme'].append('FucT')
                    d['Glycan_in'].append(add_in)
                    d['Donor'].append('GDP_Fuc')
                    d['---->'].append('---->')
                    d['Glycan_out'].append(add_out)
                    d['Residue'].append('GDP')
                while g < a: #append Galactose
                    add_in =  init_str %(0,m,a,g,0)
                    add_out = init_str %(0,m,a,g+1,0)
                    reactions_in.append(add_in)
                    reactions_out.append(add_out)
                    d['Gene'].append('GalT')
                    d['Enzyme'].append('GalT')
                    d['Glycan_in'].append(add_in)
                    d['Donor'].append('UDP_Gal')
                    d['---->'].append('---->')
                    d['Glycan_out'].append(add_out)
                    d['Residue'].append('UDP')
                    g+=1 
                    s = 0    
                    if g != 0:
                        while s < g: #append Sialic Acid
                            add_in =  init_str %(0,m,a,g,s)
                            add_out = init_str %(0,m,a,g,s+1)
                            reactions_in.append(add_in)
                            reactions_out.append(add_out)
                            d['Gene'].append('SiaT')
                            d['Enzyme'].append('SiaT')
                            d['Glycan_in'].append(add_in)
                            d['Donor'].append('CMP_Neu5AC')
                            d['---->'].append('---->')
                            d['Glycan_out'].append(add_out)
                            d['Residue'].append('CMP')
                            s+=1
    if m < 5:
        F = 1
        while a < 1: 
            add_in =  init_str %(F,m,a,0,0)
            add_out = init_str %(F,m,a+1,0,0)
            reactions_in.append(add_in)
            reactions_out.append(add_out)
            d['Gene'].append('GnTI')
            d['Enzyme'].append('GnTI')
            d['Glycan_in'].append(add_in)
            d['Donor'].append('UDP_GlcNAc')
            d['---->'].append('---->')
            d['Glycan_out'].append(add_out)
            d['Residue'].append('UDP')
            add_in =  init_str %(F,m,a+1,0,0)
            add_out = init_str %(F,m-1,a+1,0,0)
            d['Gene'].append('ManII')
            d['Enzyme'].append('ManII')
            d['Glycan_in'].append(add_in)
            d['Donor'].append('UDP_GlcNAc')
            d['---->'].append('---->')
            d['Glycan_out'].append(add_out)
            d['Residue'].append('UDP')
            a+=1
        g = 0    
        if a == 1:
            while g < 1: #append Galactose
                add_in =  init_str %(F,m,a,g,0)
                add_out = init_str %(F,m,a,g+1,0)
                reactions_in.append(add_in)
                reactions_out.append(add_out)
                d['Gene'].append('GalT')
                d['Enzyme'].append('GalT')
                d['Glycan_in'].append(add_in)
                d['Donor'].append('UDP_Gal')
                d['---->'].append('---->')
                d['Glycan_out'].append(add_out)
                d['Residue'].append('UDP')
                g+=1 
        s = 0    
        if g == 1:
            while s < 1: #append Sialic acid
                add_in =  init_str %(F,m,a,g,s)
                add_out = init_str %(F,m,a,g,s+1)
                reactions_in.append(add_in)
                reactions_out.append(add_out)
                d['Gene'].append('SiaT')
                d['Enzyme'].append('SiaT')
                d['Glycan_in'].append(add_in)
                d['Donor'].append('CMP_Neu5Ac')
                d['---->'].append('---->')
                d['Glycan_out'].append(add_out)
                d['Residue'].append('CMP')
                s+=1
    a = 0
    if m == 3: #if all manose is removed
        while a < m+1: #append GluNac m+1 indicates that 2 GlcNAc append to 1 branch of mannose
            add_in =  init_str %(F,m,a,0,0)
            add_out = init_str %(F,m,a+1,0,0)
            reactions_in.append(add_in)
            reactions_out.append(add_out)
            d['Gene'].append('GnTI')
            d['Enzyme'].append('GnTI')
            d['Glycan_in'].append(add_in)
            d['Donor'].append('UDP_GlcNAc')
            d['---->'].append('---->')
            d['Glycan_out'].append(add_out)
            d['Residue'].append('UDP')
            a+=1
            g = 0    
            if a != 0:
                while g < a: #append Galactose
                    add_in =  init_str %(F,m,a,g,0)
                    add_out = init_str %(F,m,a,g+1,0)
                    reactions_in.append(add_in)
                    reactions_out.append(add_out)
                    d['Gene'].append('GalT')
                    d['Enzyme'].append('GalT')
                    d['Glycan_in'].append(add_in)
                    d['Donor'].append('UDP_Gal')
                    d['---->'].append('---->')
                    d['Glycan_out'].append(add_out)
                    d['Residue'].append('UDP')
                    g+=1 
                    s = 0    
                    if g != 0:
                        while s < g: #append Sialic Acid
                            add_in =  init_str %(F,m,a,g,s)
                            add_out = init_str %(F,m,a,g,s+1)
                            reactions_in.append(add_in)
                            reactions_out.append(add_out)
                            d['Gene'].append('SiaT')
                            d['Enzyme'].append('SiaT')
                            d['Glycan_in'].append(add_in)
                            d['Donor'].append('CMP_Neu5AC')
                            d['---->'].append('---->')
                            d['Glycan_out'].append(add_out)
                            d['Residue'].append('CMP')
                            s+=1
            
        
df = DataFrame.from_dict(d) #dataframe from dictionary  
print(df.to_string()) #prit the dataframe of the reactions 