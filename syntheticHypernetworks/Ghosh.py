'''
Synthetic Dataset Creator, following the method in:

S.Ghosh, A.Chakraborty, "Clustering Hypergraphs for Discovery of Overlapping Communities"

@author: ak
'''
import random
import json
from math import ceil,floor
import numpy.random as nprand
import numpy as np


def Ghosh(N, L, gamma, beta, scatter=0, types=['X', 'Y', 'Z'], filepath=None,filetype='raw'):
    '''
    Synthetic Multipartite Hypergraph Dataset Creator following the method found in:
    S.Ghosh, A.Chakraborty, "Clustering Hypergraphs for Discovery of Overlapping Communities"

    allows overlapping communities
    
    Remarks:
        non overlapping communities
        does not guarrantee all nodes participate in a hyperedge
            i.e. some nodes may be unreachable
            
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
    L: int
        communities per domain/ vertex type
    gamma: float
        fraction of each partite type belonging to multiple communities
    beta: float
        hyperedge density i.e. actual_he/total_possible_he within communities
    scatter: float
        fraction of noise hyperedges
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
    
    N_multiple=int(ceil(gamma*N)) #number of nodes from each partite belonging to multiple communities
    
    nodes={t:range(1,N+1) for t in types}
    
    #partition to communities
    for t in nodes:
        random.shuffle(nodes[t])
        comms=partition(nodes[t],L)
        nodes[t]={j:comms[j] for j in range(len(comms))}
        
    comm_dict={t:{n:set([c]) for c in nodes[t] for n in nodes[t][c] } for t in types}

    #add nodes belonging to multiple communities
    for t in nodes:
        for n in random.sample(comm_dict[t].keys(),N_multiple):
            c=nprand.choice(list(set(range(L))-comm_dict[t][n]))
            comm_dict[t][n].add(c)
            nodes[t][c].append(n)
            
    M=[]
    for c in range(L):
        for t in types:
            M.append(len(nodes[t][c]))
    
    M=L*(np.mean(M)**d)
    M=int(round(beta*M)) #number of hyperedges
    noise=int(round(scatter*M)) #number of noise hyperedges
    
    
    #form hyperedges
    he=set()
    lencountcheck=0
    lenprev=0
    
    while len(he)<(M-noise):
        c=nprand.choice(range(L))
        hyperedge=[]
        for t in types:
            hyperedge.append(nprand.choice(nodes[t][c]))
        he.add(tuple(hyperedge))
        
        if lenprev==len(he):
            lencountcheck+=1
            if lencountcheck>1e4:
                #ifinite loop
                #print('BREAK FROM LOOP')
                break
        else:
            lencountcheck=0        
        lenprev=len(he)
    

    #add noise hyperedges
    while len(he)<M:
        hyperedge=[]
        for t in nodes:
            hyperedge.append(nprand.choice(nodes[t][nprand.choice(range(L))]))
        he.add(tuple(hyperedge))
            
    
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
# TEST - OK
#===============================================================================
if __name__=='__main__':
    a,b,c=Ghosh(4, 2, 0.2, 1.0, 0.2, types=['X','Y','Z', 'K'], filepath='', filetype='json')
    
    for i in a:
        print(i)
    for i in c:
        print i
        print(c[i])
    print(b)