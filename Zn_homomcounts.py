import numpy as np
import sys

data={}
dict_hash_sets={}

with open(sys.argv[1],"r") as f:
    t=f.readline()
    while(not t==""):
        t=f.readline()
        hashs = t.rstrip().replace(" ","").split("-")
        sets = f.readline(),f.readline()

        dict_hash_sets[hashs[0]]=sets[0]
        dict_hash_sets[hashs[1]]=sets[1]

        if hashs[0] in data:
            data[hashs[0]].append(hashs[1])
        else:
            data[hashs[0]]=[hashs[1]]
        if hashs[1] in data:
            data[hashs[1]].append(hashs[0])
        else:
            data[hashs[1]]=[hashs[0]]

        t=f.readline()

def follow_in_dict(the_dict,the_key,the_list):
    elems = the_dict[the_key][:]
    the_dict.pop(the_key,None)
    for x in elems:
        the_list.append(x)
        if the_dict.get(x):
            follow_in_dict(the_dict,x,the_list)

h=[]
the_keys = list(data.keys())
for k in the_keys:
    if data.get(k):
        a_list=[]
        follow_in_dict(data,k,a_list)
        a_list = np.unique(a_list)
        h.append(a_list)

print("# of homometric subsets")
val,counts = np.unique([len(x) for x in h],return_counts=True)
for x,y in zip(val,counts):
    print("{}-uples: {}".format(x,y))
