# cellCollective.py
#
# Bryan Daniels
# 2023/12/8
#
# Code for dealing with cell collective networks
#

import glob
import os
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
