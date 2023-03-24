# load_cell_collective_nets.py
#
# Bryan Daniels
# 4/12/2019
#
# Load all available networks from the Cell Collective database.
#

try:
    import datadex
except:
    print("load_cell_collective_nets WARNING: Cannot import datadex")
import numpy as np
from neet.boolean import LogicNetwork

# this version does not depend on datadex
def loadCellCollectiveNet(netdir,datadir):
    logic_file = "{}/{}/expressions.txt".format(datadir,netdir)
    external_nodes_file = "{}/{}/external.txt".format(datadir,netdir)
    return LogicNetwork.read_logic(logic_file,external_nodes_file)

# taken from modular-control.ipynb
def loadCellCollectiveNets(netindex=None,skip=['ErbB Receptor Signaling'],verbose=False):
    
    # load all cell collective networks
    # (takes about 10 seconds)
    
    masterdirectory = '../../'
    
    dexFilename = masterdirectory+'dex.db'
    dex = datadex.DataDex(dexFilename)
    
    networks = np.array( dex.select(['name','filename'],['is_biological']) )
    
    dirDict = dict(networks)
    netDict = {}
    
    if netindex is None:
        netnames = dirDict.keys()
    else:
        netnames = [ np.sort(list(dirDict.keys()))[netindex] ]
    
    for netname in netnames:
        if netname in skip:
            print("loadCellCollectiveNets: Skipping",netname)
        else:
            if verbose:
                print("loadCellCollectiveNets: Loading network",netname)
            directory = dirDict[netname]
            logic_file = masterdirectory+directory+"/expressions.txt"
            external_nodes_file = masterdirectory+directory+"/external.txt"
            net = LogicNetwork.read_logic(logic_file,external_nodes_file)
            net.metadata['name'] = netname
            netDict[netname] = net
    return netDict

def cellCollectiveNetDirectory(netindex):
    """
    Returns the datadex name and location of files for the network with index netindex.
    """
    
    masterdirectory = '../../'
    
    dexFilename = masterdirectory+'dex.db'
    dex = datadex.DataDex(dexFilename)
    
    networks = np.array( dex.select(['name','filename'],['is_biological']) )
    
    dirDict = dict(networks)

    name = np.sort(list(dirDict.keys()))[netindex]

    return name,dirDict[name]
