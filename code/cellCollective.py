# cellCollective.py
#
# Bryan Daniels
# 2023/12/8
#
# Code for dealing with cell collective networks
#

import glob
import os
import numpy as np
from neet.boolean import LogicNetwork

def load_cell_collective_network(directory):
    """
    Creates neet network from cell collective data in the
    given directory.
    """
    expressions_file = '{}/expressions.txt'.format(directory)
    external_file = '{}/external.txt'.format(directory)
    net = LogicNetwork.read_logic(expressions_file, external_file)
    return net

def load_all_cell_collective_networks(main_directory,
    skip=['ErbB_Receptor_Signaling']):
    """
    Returns dictionary of all cell collective networks
    located in subdirectories of the given main directory.
    
    skip (['ErbB_Receptor_Signaling']) : List of network names
                                         to skip
    """
    netDict = {}
    for directory in glob.glob(main_directory+"/*"):
        name = os.path.split(directory)[-1]
        if os.path.isdir(directory) and name not in skip:
            net = load_cell_collective_network(directory)
            netDict[name] = net
    return netDict

def load_cell_collective_network_from_index(main_directory,
    net_index,**kwargs):
    """
    Returns the name and neet network corresponding to the
    cell collective network with given index in the sorted
    list of names of networks located in subdirectories
    of the given main directory.
    """
    netDict = load_all_cell_collective_networks(main_directory,
        **kwargs)
    netnames = np.sort(list(netDict.keys()))
    name = netnames[net_index]
    net = netDict[name]
    return name,net
