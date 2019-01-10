#!/usr/bin/python

import pandas as pd
import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot
from string import digits
import sys
sys.tracebacklimit = 0
np.seterr(divide='ignore', invalid='ignore')

def get_atom_count(df,atom):
    atoms=[]
    for i in df[0]:
        atoms.append(" ".join([x for x in i if not x.isdigit()]))
    
    u_atom=" ".join([x for x in atom if not x.isdigit()])
    
    number=atoms.count(u_atom)
    return(number)

def get_atom_types(df):
    atoms=[]
    
    for i in df[0]:
        
        atoms.append(" ".join([x for x in i if not x.isdigit()]))
    
    atom_types=np.unique(atoms)
    
    return(atom_types)

def get_atom_groups(atom_types,df):
    #this classifies all atoms in the system into atom_groups by atom_types
    atom_groups=[[] for i in range(len(atom_types))]
    
    for i, atom in enumerate(atom_types):
        temp=[]
        for atom_index, atom_names in enumerate(df[0]):
            if atom in atom_names:
                temp.append(np.array(df.iloc[atom_index]))
        atom_groups[i].append(temp)
 
    return(atom_groups)
            #temp.append(np.array(df.iloc[atom_index]))
    
def get_group_perm(atom_groups):
    
    atom_range=list(range(0,len(atom_groups)))
    atom_group_perm=[]
    
    for p in itertools.product(atom_range, repeat=2):
        atom_group_perm.append(p)
        
    return(atom_group_perm)

def get_max_atom_distance(i,j,atom_groups):
    bond_distances=[]
    for x in itertools.product(atom_groups[i][0],atom_groups[j][0]):
        D=(x[0][1:4]-x[1][1:4])**2
        D=np.sqrt(D[0]+D[1]+D[2])
        
        bond_distances.append(D)
    #print(np.amax(bond_distances), i, j)
    return(np.amax(bond_distances))
        
        #print(x[0][0:4],x[1][0:4])
        
def get_unit_cell_diagonal(atom_range,atom_groups):
    list_of_max_bonds=[]

    for p in itertools.product(atom_range, repeat=2):
        list_of_max_bonds.append(get_max_atom_distance(p[0],p[1],atom_groups))
    
    #returns R or Diagonal/2
    return(np.amax(list_of_max_bonds))

def get_atom_density(N_atom,d):
    
    Volume=4*np.pi*(d/2)**3
    
    return(N_atom/Volume)

def get_product(i,j,atom_groups):
    for x in itertools.product(atom_groups[i][0],atom_groups[j][0]):
        print(x[0][0],x[1][0])
        #print(x[0][0:4],x[1][0:4])
        
def get_product_distance_different_atoms(i,j,atom_groups):
    for x in itertools.product(atom_groups[i][0],atom_groups[j][0]):
        X=(x[0][1]-x[1][1])**2
        Y=(x[0][2]-x[1][2])**2
        Z=(x[0][3]-x[1][3])**2
                    
        D=np.sqrt(X+Y+Z)
        
        #print(x[0][0], x[1][0], D)
        #print(x[0][0:4],x[1][0:4])


def get_RDF_histogram(i,j,atom_groups,df,d):
    
    histo_input=[]
    
    #calculate distances between atom pairs
    for x in itertools.product(atom_groups[i][0],atom_groups[j][0]):
        X=(x[0][1]-x[1][1])**2
        Y=(x[0][2]-x[1][2])**2
        Z=(x[0][3]-x[1][3])**2
                    
        D=np.sqrt(X+Y+Z)
        
        #print(x[0][0], x[1][0], D)
        histo_input.append(D)
    
    #count number of atoms in data
    #N_atom=len(histo_input)
    N_atom_1=get_atom_count(df,str(atom_groups[i][0][0][0]))
    N_atom_2=get_atom_count(df,str(atom_groups[j][0][0][0]))
    N_atom=N_atom_1+N_atom_2
    
    #get atom_density
    rho=get_atom_density(N_atom,d)
    
    #bins=int(d/0.02)
    bins=int(300)

    histo_max=d+1
    histo_min=0
    
    histo_output=np.histogram(histo_input,bins=bins,range=(histo_min,histo_max))
    
    #histo_y=histo_output[0]/(4*np.pi*histo_output[1]**2)
    histo_y=histo_output[0]/((4*np.pi*histo_output[1][:-1]**2))
    histo_y=histo_y/((rho))
    
    #normalized_histo=np.array(histo_y,histo_output[1])
    
    #print(histo_output[0],histo_output[1])
    #print(histo_output[1])
    
    ###########this is visualization of the fingerprint#############
    
    #matplotlib.pyplot.plot(histo_output[1][:-1],histo_y, 'b')
    #fig=matplotlib.pyplot.gcf()
    #fig.set_size_inches(18.5, 10.5)
    #matplotlib.pyplot.show(block=False)
    #matplotlib.pyplot.savefig('temp.png')
    
    ##########################
    
    return(histo_y)
    
    #print(len(histo_output[0]),len(histo_output[1]))

    
 
def get_fingerprint(file):
    
    fingerprint=[]
    #read xyz file in
    df=pd.read_csv(file,skiprows=[0,1], header=None, 
            delim_whitespace=True)

    #get atom_types from atom_names
    atom_types=get_atom_types(df)
    
    #get atom-atom group types
    atom_groups=get_atom_groups(atom_types,df)
    
    #get atom_atom permutation of types returns list of indexes of group_types
    atom_group_perm=get_group_perm(atom_types)
    
    #get atom_range
    atom_range=(list(range(0,len(atom_groups))))
    
    #get unit_cell diagonal
    d=get_unit_cell_diagonal(atom_range,atom_groups)
    
    for i in atom_group_perm:
        temp=get_RDF_histogram(i[0],i[1],atom_groups,df,d)
        
        temp[np.isinf(temp)] = 0
        temp[np.isnan(temp)] = 0

        fingerprint.extend(temp)

        #get_product_distance_different_atoms(i[0],i[1], atom_groups)
        #get_product_distance(i[0],i[1], atom_groups)
    
    #matplotlib.pyplot.plot(fingerprint)
    #fig=matplotlib.pyplot.gcf()
    #fig.set_size_inches(18.5, 10.5)
    #matplotlib.pyplot.show(block=False)
    
    #print(fingerprint)
    return(fingerprint)

def main():
	sys.tracebacklimit = 0
	file1=sys.argv[1]
	file2=sys.argv[2]

	histo_1=get_fingerprint(file1)
	histo_2=get_fingerprint(file2)
    
	a=np.array(np.array(histo_1[1:])+1)
	b=np.array(np.array(histo_2[1:])+1)

	signed_vector_product=sum(filter(lambda x: x != float('-inf'), a*b))
	unsigned_vector_a=sum(filter(lambda x: x != float('-inf'), abs(a)))
	unsigned_vector_b=sum(filter(lambda x: x != float('-inf'), abs(b)))
	cos_distance=signed_vector_product/(unsigned_vector_a*unsigned_vector_b)
	print(cos_distance)

if __name__=="__main__":
	main()
