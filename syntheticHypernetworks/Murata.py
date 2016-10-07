'''
Synthetic Dataset Creator, following the method in:

T.Murata, "Modularity for Heterogeneous Networks"

@author: ak
'''
from math import ceil,floor
import random
import itertools
import json

def Murata(N, R, L, p, types=['X', 'Y', 'Z'], filepath=None,filetype='raw'):
    '''
    Synthetic Multipartite Hypergraph Dataset Creator following the method found in:
    T.Murata, "Modularity for Heterogeneous Networks"
    
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
    R: float
        ratio of hyperedges per vertices
    L: int
        communities per domain/ vertex type
    p: float
        fraction(probability) of noisy nodes
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
    
    M=N*R #number of hyperedges
    
    nodes={t:range(1,N+1) for t in types}
    
    #partition to communities
    for t in nodes:
        random.shuffle(nodes[t])
        comms=partition(nodes[t],L)
        nodes[t]={j:comms[j] for j in range(len(comms))}
    
    #extract intra community hyperedges
    he=set()
    while len(he)<int(ceil(M*(1-p))):
        valid_c=random.sample(range(L),1)[0]
        he.add(tuple([random.sample(nodes[t][valid_c],1)[0] for t in types]))
    
    
    #extract noise hyperedges
    while len(he)<M:
        noise=random_permutation(range(L)*(d-1), d)
        he.add(tuple([random.sample(nodes[i][j],1)[0] for i,j in zip(types,noise)]))
        
    comm_dict={t:{n:c for c in nodes[t] for n in nodes[t][c] } for t in types}
    
    comm_array=[comm_dict[t][i] for t in types for i in sorted(comm_dict[t]) ]
    
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
                f.write(','.join([str(e) for e in comm_array]))
            
            with open(filepath+'community_matchings.raw','w') as f:
                for t in types:
                    for n in comm_dict[t]:
                        f.write(str(t)+' '+str(n)+' '+str(comm_dict[t][n])+'\n')

    return map(lambda(x):list(x),list(he)), comm_array, comm_dict
    
    
#===============================================================================
# Helper functions
#===============================================================================

def partition(lst, n):
    division = len(lst) / float(n)
    return [ lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n) ]

def random_permutation(iterable, r=None):
    "Random selection from itertools.permutations(iterable, r)"
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(random.sample(pool, r))
    
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
        
        
        
#===============================================================================
# TEST - OK
#===============================================================================
if __name__=='__main__':
    (he,comm,d)=Murata(10000, 10, 100, 0.4, types=['X','Y', 'Z','K'], filepath='.', filetype='raw')