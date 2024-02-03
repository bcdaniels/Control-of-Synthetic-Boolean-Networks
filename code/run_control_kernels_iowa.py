# run_control_kernels_iowa.py
#
# Bryan Daniels
# 2024/2/2 branched from grn-survey/scripts/modularity/run_control_kernels.py
# 2.19.2019
#
# Find control kernels for each attractor for one network in the
# Iowa database.  The network's index is set by command
# line argument.
#

import sys
from toolbox.simplePickle import save
from iowaDatabase import load_iowa_database_network_from_index
from neet.controlkernel.control_kernel_analysis import dn_analysis,ck_analysis

if len(sys.argv) != 2:
    print("Usage: python run_control_kernels_iowa.py [netIndex]")
    exit()

database_dir = "../../Data/240202/models/update_rules_models_in_literature_we_randomly_come_across/"

require_inputs = True # False # effectively True before 2020.9.17
find_control_kernels = True # False # True
find_dn_bound = False # True

# () load database network as net

netname,net = load_iowa_database_network_from_index(database_dir,
            int(sys.argv[1]) )
print("Loaded network: {}".format(netname))
print()

if find_control_kernels:
    # () run control kernel analysis
    data = ck_analysis(net)
    attractors = data["attractors"]
else:
    data = {}
    attractors = None
    
# () save output data to pickled file

data["name"] = netname
filename = "control_kernel_{}.dat".format(netname)
print("Outputting control kernel data to {}".format(filename))
print()
save(data,filename)

if find_dn_bound:
    # () run distinguishing node bound analysis using computed attractors and control kernels

    ck_list = data["control_kernels"]
    data_dn_bound = dn_analysis(net,
                                attractors=attractors,
                                possible_nodes_list=ck_list,
                                require_inputs=require_inputs)

    # () save output data to pickled file

    data.update(data_dn_bound)
    data["name"] = netname
    filename = "control_kernel_{}.dat".format(netname)
    print("Outputting distinguishing node bound data to {}".format(filename))
    print()
    save(data,filename)


# () run distinguishing node analysis

data_dn = dn_analysis(net,
                      attractors=attractors,
                      require_inputs=require_inputs)

# () save output data to pickled file

data.update(data_dn)
data["name"] = netname
filename = "control_kernel_{}.dat".format(netname)
print("Outputting distinguishing node data to {}".format(filename))
print()
save(data,filename)


