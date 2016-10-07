'''
Synthetic Dataset Creator, following the method in:

N.Neubauer, K.Obermayer, "Towards Community Detection in k-Partite, k_Uniform Hypergraph"

The "caveman" model

@author: ak
'''
import random
import json
import numpy as np
import numpy.random as nprand
import math
from __builtin__ import True

def Caveman(N, L, p, overlap=False, types=['X', 'Y', 'Z'], filepath=None,filetype='raw'):
    '''
    Synthetic Multipartite Hypergraph Dataset Creator following the method found in:
    N.Neubauer, K.Obermayer, "Towards Community Detection in k-Partite, k_Uniform Hypergraph"
    a.k.a. the "caveman" model
    
    Remarks:
        allows overlapping communities
        allows different number of communities for each partite
        
    In case of json files:
    (a) hyperedges.json file: {"hyperedges": [[<idXtype> , <idYtype>, ..., <idKtype>],
                                                ...
                                              [<idXtype> , <idYtype>, ..., <idKtype>]}
    (b)community_array.json file: {"comms": [<commId1>, ..., <commIdN>]} 
                    , communities of nodes ordered by partite, by id
    (c)community_matchings.json file {<partiteId>:{<nodeId>:<communityId>,
                                                    ...
                                                  },
                                      ...,
                                      <partiteKId>:{<nodeId>:<communityId>,
                                                    ...
                                                  }
                                    }
    
    In case of raw files:
    (a) hyperedges.raw file: <idXtype> <idYtype> ... <idKtype>
                             ...
                             <idXtype> <idYtype> ... <idKtype>
                space delimited node ids,
                one hyperedge per line
    (b)community_array.raw file: <commId1>, ..., <commIdN>
                communities of nodes ordered by partite, by id, separated by comma (,)
    (c)community_matchings.raw file: <partiteId> <nodeId> <communityId>
                                     <partiteId> <nodeId> <communityId>
                                     ...
                                     <partiteId> <nodeId> <communityId>
                one line per node,
                the partite it belongs to, the node id, the community id it belongs to

    Parameters            
    ----------
    N: int
        number of vertices in each domain
    L: int or list
        communities per domain/ vertex type
    p: float
        probability 1-(1-p)^(1/(k-1)) of hyperedge within corresponding communities
    overlap: bool, optional
        if True allow overlapping communities
        default=False
    types: list of str, optional
        domain names/vertex types in multipartite network, default=['X', 'Y', 'Z']
        types also define the k in the k-partite network
    filepath: str, optional
        if specified saves the dataset to given path
        produces 3 files : (a) community_array (b) hyperedges (c) vertex_comm_dictionary
    filetype: str, optional
        defines format of produced files
         'raw', 'json'
         
    Returns
    -------
    he: list of lists
        list of hyperedges, each hyperedge a list of vertices
    comm_array: list
        community array, in order of <types> and in sorted vertex id order
    comm_dict: dictionary
        dictionary of vertex-communityId pairs
    '''
    d=len(types) #number of domains/partites
    
    if isinstance(L,int):
        L=[L]*d
        
    nodes={t:range(1,N+1) for t in types}

    
    #partition to communities
    for t in nodes:
        random.shuffle(nodes[t])
        comms=partition(nodes[t],L[types.index(t)])
        nodes[t]={j:comms[j] for j in range(len(comms))}
    
    for t in types:
        print(t)
        print(nodes[t])
        
    nodes_tupled=[]
    for t in types:
        for c in nodes[t]:
            nodes_tupled.extend([(t, n, c) for n in nodes[t][c]])
    random.shuffle(nodes_tupled)

    comm_dict={t:{n:set([c]) for c in nodes[t] for n in nodes[t][c] } for t in types}
        
    print(nodes_tupled)
    #build hyperedges
    he=set()
    for nt,nId,nc in nodes_tupled:
        print(str(nt)+' '+str(nId)+' '+str(nc))
        hyperedge=[]
        he_existed=False
        while (he_existed==True) or (he_existed==False and len(hyperedge)==0):
            
            #check extreme condition, all possible combinations exist
            if len(he)==N**d:
                break
            
            for t in types:
                if t==nt:
                    hyperedge.append(nId)
                else:
                    if overlap:
                        probs=np.array(map(lambda x: 1.0-math.pow(1.0-p, 1.0/(d-1.0)) if x==(nc%L[types.index(t)]) else math.pow(1.0-p, 1.0/(d-1.0)), 
                                               range(L[types.index(t)])))
                        if probs.sum()!=1.0:
                            probs/=probs.sum()
                            
                        #===========================================================
                        # print('-----------------------------')
                        # print((nt,nId,nc))
                        # print(t)
                        # print(range(L[types.index(t)]))
                        # print(probs)
                        # print(nc%L[types.index(t)])
                        # print('-----------------------------')
                        #===========================================================
                            
                        c=nprand.choice(range(L[types.index(t)]),p=probs)
                    else:
                        probs=np.array(map(lambda x: 1.0-math.pow(1.0-p, 1.0/(d-1.0)) if x==nc else math.pow(1.0-p, 1.0/(d-1.0)), 
                                               range(L[types.index(t)])))
                        if probs.sum()!=1.0:
                            probs/=probs.sum()
                        c=nprand.choice(range(L[types.index(t)]),p=probs)
                    
                    mn=nprand.choice(nodes[t][c])
                    
                    if overlap and c==(nc%L[types.index(t)]):
                        #node belongs to overlapping communities
                        comm_dict[nt][nId].add(c)
                    hyperedge.append(mn)
            
            if (tuple(hyperedge) in he):
                he_existed=True
                hyperedge=[]
            else:
                he.add(tuple(hyperedge))
                break
        
    for t in types:
        for n in comm_dict[t]:
            comm_dict[t][n]=list(comm_dict[t][n])
    
    comm_array=[list(comm_dict[t][i]) for t in types for i in sorted(comm_dict[t]) ]
    
    if filepath:
        if filetype=='json':
            json_writer({'hyperedges':map(lambda(x):list(x),list(he))}, filepath+'hyperedges.json')
            json_writer({'comms':comm_array}, filepath+'community_array.json')
            json_writer(comm_dict, filepath+'community_matchings.json')
        else:
            with open(filepath+'hyperedges.raw','w') as f:
                for e in he:
                    f.write(' '.join([str(i) for i in e])+'\n')
            
            with open(filepath+'community_array.raw','w') as f:
                f.write('\n'.join([list_to_raw_string(e) for e in comm_array]))
            
            with open(filepath+'community_matchings.raw','w') as f:
                for t in types:
                    for n in comm_dict[t]:
                        f.write(str(t)+' '+str(n)+' '+list_to_raw_string(comm_dict[t][n])+'\n')

    return map(lambda(x):list(x),list(he)), comm_array, comm_dict
        
        
#===============================================================================
# Helper functions
#===============================================================================

def partition(lst, n):
    division = len(lst) / float(n)
    return [ lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n) ] 
        
def json_writer(obj, filename):
    """write obj data (assume well formatted) to json file
    
    Parameters
    ----------
    obj: object to dump to json
    filename: string
        filename of json file (with extension .json)
    
    """
    with open(filename, 'w') as f:
        json.dump(obj, f)
        
def list_to_raw_string(l):
    s=str(l[0])
    for i in l[1:]:
        s+=','+str(i)
    return s
#===============================================================================
# TEST
#===============================================================================
if __name__=='__main__':
    a,b,c=Caveman(2, [1,1,1,1], 0.95, overlap=False, types=['X','Y','Z','K'], filepath='', filetype='json')
    
    for i in a:
        print(i)
    for i in c:
        print i
        print(c[i])
    print(b)