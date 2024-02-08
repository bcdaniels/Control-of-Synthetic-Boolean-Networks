# run_basin_data_iowa.py
#
# Bryan Daniels
# 2024/2/8
#
# compute is_isolated_list for Iowa database networks (not bothering for now with using sampled basins)

from load_control_kernel_data import loadDataExact
from iowaDatabase import load_all_iowa_database_networks
from preimages import isolated_list
from toolbox.simplePickle import save

def computeBasinData(dataDict,netDict):
    basinData = {}
    for name in dataDict:
        print("run_basin_data_iowa: Analyzing {}...".format(name))
        CKdata = dataDict[name]
        net = netDict[name]
        is_isolated_list = isolated_list(net,CKdata['attractors'])
        basinData[name] = {'is_isolated_list':is_isolated_list,
                           'attractors':CKdata['attractors'],
                                  }
    return basinData
    
if __name__=='__main__':
    iowaNetDir = '/Users/bdaniel6/ASUDropbox/Research/control-continuous/Data/240202/models/update_rules_models_in_literature_we_randomly_come_across/'
    iowaDataDir = '/Users/bdaniel6/ASUDropbox/Research/control-continuous/Data/240203-control-kernels/'
    netDict = load_all_iowa_database_networks(iowaNetDir)
    dataDictIowa,dfIowa = loadDataExact(iowaDataDir)
    basinDataIowa = computeBasinData(dataDictIowa,netDict)
    save(basinDataIowa,'240208_basinDataIowa.pkl')
