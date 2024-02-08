# preimages.py
#
# Bryan Daniels, Enrico Borriello
#
# Functions for computing the preimages of any given state
# in a Neet network.
#

import pandas as pd
from neet import UniformSpace
import tqdm

def activating_states(net,node_index,activate=True):
    """
    Returns the states that activate or deactivate a given node in a Neet network.
    
    Returns: (neighbors_tuple,conditions),
    where neighbors_tuple is a sorted tuple of node indices on which the given node depends,
    and conditions is a tuple of binary states that activate the given node
    (or deactivate the given node if activate=False)
    """
    neighbors_list = tuple(sorted(net.neighbors_in(node_index)))
    conditions = []
    default_state = [ 0 for i in range(net.size) ] # states to which other non-input nodes are set
    input_space = UniformSpace(len(neighbors_list),2) # create binary input space
    for reduced_state in input_space:
        
        # set starting state to have given values for input nodes
        full_state = default_state.copy()
        for i,neighbor in enumerate(neighbors_list):
            full_state[neighbor] = reduced_state[i]
        
        # check whether this state is activating for node_index
        if net.update(full_state)[node_index] == activate:
            conditions.append(tuple(reduced_state))
            
    return neighbors_list,tuple(conditions)

# from https://github.com/EnricoBorriello/Boolean-Backwards-Reachability/blob/main/backreach.py
def conditions_product(df1,df2):
    '''
    Activating conditions for functions 1 AND 2 given the input truth tables
    for function 1 and function 2. The inputs need to be dataframes with column names
    corresponding to the the labes of the nodes.
    '''

    # Firts check if the two function share input nodes

    shared_columns = list(set(df1.columns) & set(df2.columns))

    if len(shared_columns) == 0:
        #if they don't, take all possible products of conditions
        merged_df = pd.merge(df1, df2, how = 'cross')
    else:
        #if they do, just keep the non conflicting ones
        merged_df = pd.merge(df1, df2, on = shared_columns)

    merged_df = merged_df.sort_index(axis=1)

    return merged_df

def conditions_product_list(df_list):
    """
    Return the product of a list of condition dataframes.
    """
    df = df_list[0]
    for dfi in df_list[1:]:
        df = conditions_product(df,dfi)
    return df

def activating_conditions_df(net,node_index,node_state):
    nodes,conditions = activating_states(net,node_index,activate=node_state)
    return pd.DataFrame(conditions,columns=nodes)

def preimages(net,state,communities=None):
    """
    Given a state of a Neet network, returns a list of all primages of that state
    (all states that produce the given state after one timestep)
    """
    if communities == None:
        communities = [range(net.size)]
    
    df_list = []
    for community in communities:
        community = list(community)
        df = activating_conditions_df(net,community[0],state[community[0]])
        for i in range(1,len(community)):
            dfi = activating_conditions_df(net,community[i],state[community[i]])
            df = conditions_product(df,dfi)
        df_list.append(df)
    #return df_list
    return conditions_product_list(df_list)

def isolated_list(net,attractors,basin_samples=None):
    """
    Returns a list of Boolean values of length number of attractors
    corresponding to whether each attractor is "isolated"
    (has basin of size 1).
    
    basin_samples (None)        : Optionally give list of basin
                                  samples to avoid computing
                                  preimages of attractors that
                                  we already know are not
                                  isolated.
    """
    if basin_samples is None:
        basin_samples = [ 0 for a in attractors ]
    assert(len(attractors)==len(basin_samples))
    
    is_isolated_list = []
    for i,att in enumerate(tqdm.tqdm(attractors)):
        if len(att) > 1 or basin_samples[i] > 1:
            is_isolated_list.append(False)
        else:
            decoded_att = net.decode(att[0])
            is_isolated_list.append(
                len(preimages(net,decoded_att))==1 )
    return is_isolated_list
