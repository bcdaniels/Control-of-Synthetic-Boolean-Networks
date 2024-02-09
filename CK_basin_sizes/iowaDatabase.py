# iowaDatabase.py
#
# Bryan Daniels
# 2024/2/2 branched from cellCollective.py
# 2023/12/8
#
# Code for dealing with networks in the Iowa database
#

import glob
import os
import numpy as np
from neet.boolean import LogicNetwork

def load_iowa_database_network(expressions_file,external_file):
    """
    Creates neet network from iowa database data in the
    given directory.
    """
    net = LogicNetwork.read_logic(expressions_file,
        external_file,parentheses_in_names=False)
    return net

def load_all_iowa_database_networks(directory,
    skip=['27765040_tabular','26603105_tabular',
    '33780439','23658556_model_10','23935937']):
    """
    Returns dictionary of all iowa database networks
    located in the given directory.
    
    skip (['27765040_tabular',   : List of network names to skip
           '26603105_tabular',   : (uninterpretable files)
           '33780439',           : |
           '23658556_model_10']) : /
           '23935937'])          : <- too many inputs to one gene
    """
    netDict = {}
    for filename in glob.glob(directory+"/*.txt"):
        name = os.path.split(filename)[-1][:-4]
        if os.path.isdir(directory) \
            and name not in skip    \
            and not name.endswith('external'):
            external_filename = filename.replace('.txt',
                                                 '_external.txt')
            net = load_iowa_database_network(filename,
                                             external_filename)
            netDict[name] = net
    return netDict
    
def load_iowa_database_network_from_index(directory,
    net_index,**kwargs):
    """
    Returns the name and neet network corresponding to the
    iowa database network with given index in the sorted
    list of names of networks located in the given
    directory.
    """
    netDict = load_all_iowa_database_networks(directory,
        **kwargs)
    netnames = np.sort(list(netDict.keys()))
    name = netnames[net_index]
    net = netDict[name]
    return name,net

def find_external_nodes(expressions_file):
    """
    Load iowa database network file and find the names of
    all external nodes.  We assume any name that does not
    appear on the left hand side of an equation is external.
    """
    
    # Uses code from neet.boolean.LogicNetwork.read_logic
    
    names = []
    expressions = []
    with open(expressions_file) as eq_file:
        for eq in eq_file:
            # ensure parentheses are separated from names
            eq = eq.replace('(',' ( ').replace(')',' ) ').strip()
            if len(eq) > 0: # skip blank lines
                name, expr = eq.split('=')
                names.append(name.strip())
                expressions.append(expr.strip())
                
    ops = {'AND', 'OR', 'NOT'}
    
    extras = []
    for expr in expressions:
        sub_nodes = []
        conditions = set()

        expr_split = expr.split()
        for i, item in enumerate(expr_split):
            if item not in ops and item not in '()':
                if item not in names:
                    extras.append(item)
        
    return set(extras)

def write_external_nodes_file(expressions_file):
    external_nodes = find_external_nodes(expressions_file)
    external_file = expressions_file[:-4] + '_external.txt'
    with open(external_file,'w') as outfile:
        outfile.write('\n'.join(str(n) for n in external_nodes))

def write_all_external_nodes_files(directory,
    skip=['27765040_tabular','26603105_tabular']):
    """
    Needs to be run once to write external nodes files
    for all networks.
    """
    for filename in glob.glob(directory+"/*.txt"):
        name = os.path.split(filename)[-1][:-4]
        if os.path.isdir(directory) and not name.endswith('external') and name not in skip:
            write_external_nodes_file(filename)

