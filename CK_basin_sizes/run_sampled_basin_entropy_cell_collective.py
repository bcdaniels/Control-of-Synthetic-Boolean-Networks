# run_sampled_basin_entropy_cell_collective.py
#
# Bryan Daniels
# 2023/12/8
#
# Calculate sampled basin entropies for cell collective networks.
#
# Runs computations for a single network based on a command
# line argument specifying the requested network index.
#

import csv
import sys
from datetime import datetime
from toolbox.simplePickle import save
from cellCollective import load_cell_collective_network_from_index
import neet.controlkernel.modularity as md
from InfEst.entropyEstimates import meanAndStdevEntropyNem

# set parameters
num_samples = 10000
cell_collective_directory = '../../Data/Cell Collective/'

# set network to run from command line argument
if len(sys.argv) != 2:
    print("Usage: python run_sampled_basin_entropy_cell_collective.py [net_index]")
    exit()
net_index = int(sys.argv[1])
name,net = load_cell_collective_network_from_index(
            cell_collective_directory,net_index)
print("run_sampled_basin_entropy_cell_collective: "\
      "Running analysis for network {}...".format(name))
      
# set output filenames
csv_filename = 'basin_entropy_data_{}.csv'.format(name)
data_filename = 'basin_entropy_data_{}.pkl'.format(name)

with open(csv_filename,mode='w',newline='',buffering=1) as file:
    
    # set up output file
    writer = csv.writer(file)
    # write header to output file
    writer.writerow(['name',
                     'network_size',
                     'num_states',
                     'num_samples',
                     'num_attractors',
                     'basin_entropy_NSB',
                     'std_basin_entropy_NSB',
                     'elapsed_time'])
                     
    start_time = datetime.now()
                    
    # do basin entropy analysis
    attractors,freqs = \
        md.sampled_basin_counts(net,num_samples=num_samples)
        
    basin_entropy_NSB,std_basin_entropy_NSB = \
        meanAndStdevEntropyNem(freqs)
                        
    end_time = datetime.now()
    elapsed_time = end_time - start_time
                        
    network_size = net.size
    num_states = 2**network_size
    num_attractors = len(attractors)

    # write data to output csv file
    writer.writerow(
        [name,
         network_size,
         num_states,
         num_samples,
         num_attractors,
         basin_entropy_NSB,
         std_basin_entropy_NSB,
         elapsed_time])
                    
    # write more detailed data to output pickle file
    dataDict = {'name':name,
                'network_size':network_size,
                'num_states':num_states,
                'attractors':attractors,
                'num_attractors':num_attractors,
                'num_samples':num_samples,
                'sampled_basin_freqs':freqs,
                'basin_entropy_NSB':basin_entropy_NSB,
                'std_basin_entropy_NSB':std_basin_entropy_NSB,
                'elapsed_time':elapsed_time,
                }
    save(dataDict,data_filename)
