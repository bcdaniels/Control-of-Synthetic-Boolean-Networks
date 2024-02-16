# load_control_kernel_data.py
#
# Bryan Daniels
# 2024/1/5 copied from grn-survey/scripts/modularity
# 4.11.2019
#
# Load data created by run_control_kernel.py.
#

import pandas as pd
import numpy as np
import glob
from toolbox.simplePickle import load
#import datadex # for cellCollectiveNameFromID

def cellCollectiveNameFromID(masterDatadexDirectory='../../'):
    dexFilename = masterDatadexDirectory+'dex.db'
    dex = datadex.DataDex(dexFilename)
    networkNames = np.array( dex.select(['name',],['is_biological']) )
    return dict(enumerate(np.sort(networkNames)))

def cellCollectiveSize(masterDatadexDirectory='../../'):
    dexFilename = masterDatadexDirectory+'dex.db'
    dex = datadex.DataDex(dexFilename)
    networkNamesAndsizes = np.array( dex.select(['name','node_count'],['is_biological']) )
    return pd.Series(dict(networkNamesAndsizes))

def _error_file_summary(filename):
    with open(filename) as f:
        for line in f:
            pass
        last_line = line[:25].strip() # (don't include line ending character)
    return last_line

def loadErrorData(dir='.',verbose=True,includeSize=True):
    """
    Get summary information about errors and associate with network names.
    """
    errDict = {}
    nameFromID = cellCollectiveNameFromID()
    for id,name in nameFromID.items():
        problem = False
        error_files = glob.glob("{}/*_{}.err".format(dir,id))
        if len(error_files) == 1:
            try:
                summary = _error_file_summary(error_files[0])
            except:
                problem = True
        else:
            problem = True
        if problem:
            if verbose:
                print("loadErrorData: Could not load error data for {} {}".format(id,name))
        else:
            errDict[name] = id,summary
    
    df = pd.DataFrame.from_dict(errDict,orient='index',columns=['id','error'])
    if includeSize:
        df['size'] = cellCollectiveSize()
    return df

def loadErrorDataRandom(dir='.'):
    errDict = {}
    for error_filename in glob.glob("{}/slurm.rbn_control_kernels.*.err".format(dir)):
        try:
            summary = _error_file_summary(error_filename)
            id = error_filename.split('.')[-2].split('_')[-1]
            errDict[id] = summary
        except:
            print("loadErrorDataRandom: Could not process {}".format(filename))
    
    df = pd.DataFrame.from_dict(errDict,orient='index',columns=['error'])
    return df
  
def loadDataModules(dir='.'):
    dataDict = {}
    dataDictForFrame = {}
    for filename in glob.glob(dir+'/modules_*.dat'):
        try:
            d = load(filename)
            success = True
        except:
            print("loadDataModules: Error loading file {}".format(filename))
            success = False
            
        if success:
            # attempting to cure problems with loading python 2 data
            for key in list(d.keys()):
                if type(key) == bytes:
                    d[key.decode()] = d[key]
            
            dataDict[d['name']] = d
            
            modules = d['modules']
            
            ddata = {'number of modules': len(modules),
                     'average module size': np.mean([ len(m) for m in modules ]),
                     'max module size': np.max([ len(m) for m in modules ]),
                     'min module size': np.min([ len(m) for m in modules ]),
            }
            
            dataDictForFrame[d['name']] = ddata
            
    df = pd.DataFrame.from_dict(dataDictForFrame,'index')

    return dataDict,df

def loadDataBasins(dir='.'):
    dataDict = {}
    dataDictForFrame = {}
    for filename in glob.glob(dir+'/basins_*.dat'):
        try:
            d = load(filename)
            success = True
        except:
            print("loadDataBasins: Error loading file {}".format(filename))
            success = False
            
        if success:
            # attempting to cure problems with loading python 2 data
            for key in list(d.keys()):
                if type(key) == bytes:
                    d[key.decode()] = d[key]
            
            dataDict[d['name']] = d
            
            basin_sizes = d['basin_sizes']
            
            if len(basin_sizes) > 1:
                maxgap = max(np.sort(basin_sizes)[1:] - np.sort(basin_sizes)[:-1])
            else:
                maxgap = 0
            
            ddata = {'number of attractors': len(basin_sizes),
                     'average basin size': np.mean(basin_sizes),
                     'std basin size': np.std(basin_sizes),
                     'max basin size': np.max(basin_sizes),
                     'min basin size': np.min(basin_sizes),
                     'basin entropy': d['basin_entropy'],
                     'number of size one basins': np.sum(basin_sizes == 1),
                     'max gap': maxgap,
                     'basin time': d['basin_time_minutes'],
            }
            
            dataDictForFrame[d['name']] = ddata
            
    df = pd.DataFrame.from_dict(dataDictForFrame,'index')

    return dataDict,df

def loadDataExact(dir='.',require_ck=False):
    dataDict = {}
    dataForFrame = []
    for filename in glob.glob(dir+'/control_kernel_*.dat'):
        if filename.find('split') == -1: # we don't want to include "split" data
            try:
                d = load(filename)
                success = True
            except:
                print("loadDataExact: Error loading file {}".format(filename))
                success = False
        else:
            success = False
            
        if success:
            # attempting to cure problems with loading python 2 data
            for key in list(d.keys()):
                if type(key) == bytes:
                    d[key.decode()] = d[key]
        
            if (not require_ck) or ('control_kernels' in d):
                dataDict[d['name']] = d
            
            ddata = dataFrameExact(d,require_ck=require_ck)
            
            if len(ddata) > 0:
                dataForFrame.append(ddata)
                
    df = pd.DataFrame.from_records(dataForFrame)
    df.set_index('name',inplace=True)
    return dataDict,df

        

def dataFrameExact(d,require_ck=False):
    ddata = {'name':d['name']}
            
    if 'encoded_input' in d:
        ddata['encoded input'] = d['encoded_input']
        
    if 'threshold parameter' in d:
        ddata['threshold'] = d['threshold parameter']
        
    if 'mean degree' in d:
        ddata['mean degree'] = d['mean degree']
        
    if 'inhibitory probability' in d:
        ddata['inhibitory probability'] = d['inhibitory probability']
        
    if 'network seed' in d:
        ddata['network seed'] = d['network seed']
            
    if 'attractors' in d:
        fixedPointIndices = [ i for i in range(len(d['attractors'])) \
                                  if len(d['attractors'][i]) == 1 ]
        attractorsFilteredFixedPoint = [ d['attractors'][i] \
                                         for i in fixedPointIndices ]
        
        ddata.update({'size': d['size'],
                 'has limit cycles': d.get('has_limit_cycles',None),
                 'number of attractors': len(d['attractors']),
                 'number of fixed point attractors': len(attractorsFilteredFixedPoint),
                 'is cell collective network': not d['name'].startswith('random'),
                 })
    
    if require_ck and ('control_kernels' not in d):
        return {}
    
    if 'control_kernels' in d:
        # get control kernel data
        
        # note: we restrict the analysis to attractors that have
        # control kernels
        indicesWithControlKernel = [ i for i in range(len(d['attractors'])) \
                                     if d['control_kernel_sizes'][i] is not None ]
        sizesFiltered = [ d['control_kernel_sizes'][i] \
                          for i in indicesWithControlKernel ]
        attractorsFiltered = [ d['attractors'][i] for i in indicesWithControlKernel ]
        
        # separately restrict to fixed point attractors (no cycles)
        sizesFilteredFixedPoint = [ d['control_kernel_sizes'][i] \
                                    for i in fixedPointIndices ]
        
        # check if control kernels for all attractors are identical
        ck_identical = np.all([ ck == d['control_kernels'][0] for ck in d['control_kernels'] ])
        
        ddata.update({
                 'mean control kernel size': np.mean(sizesFiltered),
                 'std control kernel size': np.std(sizesFiltered),
                 'mean fixed point control kernel size': np.mean(sizesFilteredFixedPoint),
                 'std fixed point control kernel size': np.std(sizesFilteredFixedPoint),
                 'has simple control entropies': d.get('simple_control_entropies',None),
                 'number of attractors with control kernel': len(attractorsFiltered),
                 'all control kernels are identical': ck_identical,
                 'control kernel time': d.get('control_kernel_time_minutes',None),
                })
                
        if 'delta_control_nodes' in d:
            ddata['number of modules'] = len(d['delta_control_nodes'])
    # end control kernel data
    
    if 'distinguishing_nodes' in d:
        # get distinguishing node data
        
        # note: we restrict the analysis to attractors that have CONTROL KERNELS
        # (there are cases with distinguishing nodes but no control kernel, which
        #  we now filter out)
        
        if 'control_kernels' not in d:
            print("WARNING: 2020.9.17 For now, taking mean over attractors with dist. nodes instead of those with control kernels")
            print("WARNING: 2020.9.17 Network name = {}".format(d['name']))
            # THE FOLLOWING LINE SHOULD BE REMOVED IN THE FINAL ANALYSIS
            indicesWithControlKernel = [ i for i in range(len(d['attractors'])) \
                                         if d['distinguishing_nodes_sizes'][i] is not None ]
        
        dnSizesFiltered = [ d['distinguishing_nodes_sizes'][i] \
                          for i in indicesWithControlKernel ]
        attractorsFilteredDN = [ d['attractors'][i] for i in indicesWithControlKernel ]
        
        ddata.update({
            'number of attractors with distinguishing nodes': len(attractorsFilteredDN),
            'mean number of distinguishing nodes': np.mean(dnSizesFiltered),
            'std number of distinguishing nodes': np.std(dnSizesFiltered),
            'distinguishing nodes time': d.get('distinguishing_nodes_time_minutes',None),
            })
            
    if 'distinguishing_nodes_with_inputs' in d:
        # get distinguishing node data
        
        # note: we restrict the analysis to attractors that have CONTROL KERNELS
        # (there are cases with distinguishing nodes but no control kernel, which
        #  we now filter out)
        
        if 'control_kernels' not in d:
            print("WARNING: 2020.9.17 For now, taking mean over attractors with dist. nodes instead of those with control kernels")
            print("WARNING: 2020.9.17 Network name = {}".format(d['name']))
            # THE FOLLOWING LINE SHOULD BE REMOVED IN THE FINAL ANALYSIS
            indicesWithControlKernel = [ i for i in range(len(d['attractors'])) \
                                         if d['distinguishing_nodes_with_inputs_sizes'][i] is not None ]
        
        dnSizesFiltered = [ d['distinguishing_nodes_with_inputs_sizes'][i] \
                          for i in indicesWithControlKernel ]
        attractorsFilteredDN = [ d['attractors'][i] for i in indicesWithControlKernel ]
        
        ddata.update({
            'number of attractors with distinguishing nodes with inputs': len(attractorsFilteredDN),
            'mean number of distinguishing nodes with inputs': np.mean(dnSizesFiltered),
            'std number of distinguishing nodes with inputs': np.std(dnSizesFiltered),
            'distinguishing nodes time': d.get('distinguishing_nodes_time_minutes',None),
            })
    
    if 'distinguishing_nodes_bound' in d:
        # get distinguishing node bound data
        
        # note: we restrict the analysis to attractors that have CONTROL KERNELS
        # (there are cases with distinguishing nodes but no control kernel, which
        #  we now filter out)
        #indicesWithDistNodesB = [ i for i in range(len(d['attractors'])) \
        #                             if d['distinguishing_nodes_bound_sizes'][i] is not None ]
        dnBSizesFiltered = [ d['distinguishing_nodes_bound_sizes'][i] \
                          for i in indicesWithControlKernel ]
        attractorsFilteredDNB = [ d['attractors'][i] for i in indicesWithControlKernel ]
        
        ddata.update({
            'number of attractors with distinguishing nodes bound': len(attractorsFilteredDNB),
            'mean number of distinguishing nodes bound': np.mean(dnBSizesFiltered),
            'std number of distinguishing nodes bound': np.std(dnBSizesFiltered),
            'distinguishing nodes bound time': d.get('distinguishing_nodes_bound_time_minutes',None),
            })
            
    return ddata

def loadDataSampled(dir='.'):
    dataDict = {}
    dataForFrame = [] #dataDictForFrame = {}
    for filename in glob.glob(dir+'/sampled_control_kernel_*.dat'):
        try:
            d = load(filename)
            success = True
        except:
            print("loadDataSampled: Error loading file {}".format(filename))
            success = False
            
        if success:
            d = load(filename)
            dataDict[d['name']] = d
            
            ddata = dataFrameSampled(d)
            
            if len(ddata) > 0:
                dataForFrame.append(ddata) # dataDictForFrame[d['name']] = ddata
            
    # end loop over networks

    df = pd.DataFrame.from_records(dataForFrame) #,'index') # from_dict
    df.set_index('name',inplace=True)

    return dataDict,df


def dataFrameSampled(d):
    ddata = {'name':d['name']}
    
    # 8.12.2019 load distinguishing node data
            
    #dnSizes = d.get( 'sampled_distinguishing_nodes_sizes', [] )
    #dnSizesFiltered = [ dnSizes[i] for i in indicesWithControlKernel ]
    
    # note: we restrict the dist. node analysis to attractors that have
    # static dist. nodes (this can in some cases be larger than the
    # number of attractors with static control kernels)
    
    if 'encoded_input' in d:
        ddata['encoded input'] = d['encoded_input']
    
    if 'sampled_distinguishing_nodes_with_inputs' in d:
    
        dnSizesWInputs = d.get( 'sampled_distinguishing_nodes_with_inputs_sizes', [] )
        indicesWithDistNodesWInputs = [ i for i in range(len(dnSizesWInputs)) \
                                        if dnSizesWInputs[i] is not None ]
        dnSizesWInputsFiltered = [ dnSizesWInputs[i] for i in indicesWithDistNodesWInputs ]
        numAttractors = len(dnSizesWInputs)
        numAttractorsFiltered = len(dnSizesWInputsFiltered)
        
        ddata.update( {'size': d['size'],
             'number of attractors': numAttractors,
             'number of attractors with distinguishing nodes': numAttractorsFiltered,
             'mean number of distinguishing nodes with inputs': np.mean(dnSizesWInputsFiltered),
             'std number of distinguishing nodes with inputs': np.std(dnSizesWInputsFiltered),
             'distinguishing node time': d.get('sampled_distinguishing_nodes_time_minutes',None),
             'is cell collective network': not d['name'].startswith('random'),
            } )
    elif 'sampled_distinguishing_nodes_sizes' in d: # inputs not included
        dnSizes = d.get( 'sampled_distinguishing_nodes_sizes', [] )
        indicesWithDistNodes = [ i for i in range(len(dnSizes)) \
                                     if dnSizes[i] is not None ]
        dnSizesFiltered = [ dnSizes[i] for i in indicesWithDistNodes ]
        numAttractors = len(dnSizes)
        numAttractorsFiltered = len(dnSizesFiltered)
    
        ddata.update({'size': d['size'],
             'number of attractors': numAttractors,
             'number of attractors with distinguishing nodes': numAttractorsFiltered,
             'mean number of distinguishing nodes': np.mean(dnSizesFiltered),
             'std number of distinguishing nodes': np.std(dnSizesFiltered),
             'distinguishing node time': d.get('sampled_distinguishing_nodes_time_minutes',None),
             'is cell collective network': not d['name'].startswith('random'),
            })
    else: # no distinguishing node data
        numAttractors = None
        ddata.update({'size': d['size'],
            'is cell collective network': not d['name'].startswith('random'),
            })
    # end distinguishing node data
    
    # record control kernel data if we have it
    if 'sampled_control_kernel_sizes' in d:
        # note: we restrict the analysis to attractors that have control kernels
        indicesWithControlKernel = [ i for i in range(len(d['sampled_attractors'])) \
                                     if d['sampled_control_kernel_sizes'][i] is not None ]
        sizesFiltered = [ d['sampled_control_kernel_sizes'][i] for i in indicesWithControlKernel ]
        if 'sampled_control_kernel_iterative_sizes' in d:
            sizesIterativeFiltered = [ d['sampled_control_kernel_iterative_sizes'][i] for i in indicesWithControlKernel ]
            iterativeRoundsFiltered = [ len(d['sampled_iterative_rounds_list'][i]) for i in indicesWithControlKernel ]
        else:
            sizesIterativeFiltered = [ np.nan for i in indicesWithControlKernel ]
            iterativeRoundsFiltered = [ np.nan for i in indicesWithControlKernel ]
        attractorsFiltered = [ d['sampled_attractors'][i] for i in indicesWithControlKernel ]
        
        if numAttractors:
            # check for consistency with distinguishing node analysis
            assert(numAttractors == len(d['sampled_attractors']))
        
        # separately restrict to fixed point attractors (no cycles)
        indicesWithFixedPoint = [ i for i in range(len(d['sampled_attractors'])) \
                                  if len(d['sampled_attractors'][i]) == 1 ]
        sizesFilteredFixedPoint = [ d['sampled_control_kernel_sizes'][i] \
                                    for i in indicesWithFixedPoint ]
        if 'sampled_control_kernel_iterative_sizes' in d:
            sizesIterativeFilteredFixedPoint = [ \
                                        d['sampled_control_kernel_iterative_sizes'][i] \
                                        for i in indicesWithFixedPoint ]
            iterativeRoundsFilteredFixedPoint = [ len(d['sampled_iterative_rounds_list'][i]) \
                                        for i in indicesWithFixedPoint ]
        else:
            sizesIterativeFilteredFixedPoint = [ np.nan for i in indicesWithFixedPoint ]
            iterativeRoundsFilteredFixedPoint = [ np.nan for i in indicesWithFixedPoint ]
        attractorsFilteredFixedPoint = [ d['sampled_attractors'][i] \
                                         for i in indicesWithFixedPoint ]
        
        # check if control kernels for all attractors are identical
        ck_identical = np.all([ ck == d['sampled_control_kernels'][0] for ck in d['sampled_control_kernels'] ])
    
        ddataCK = {
             'size': d['size'],
             'mean control kernel size': np.mean(sizesFiltered),
             'std control kernel size': np.std(sizesFiltered),
             'mean iterative control kernel size': np.mean(sizesIterativeFiltered),
             'std iterative control kernel size': np.std(sizesIterativeFiltered),
             'mean fixed point control kernel size': np.mean(sizesFilteredFixedPoint),
             'std fixed point control kernel size': np.std(sizesFilteredFixedPoint),
             'mean iterative fixed point control kernel size': np.mean(sizesIterativeFilteredFixedPoint),
             'std iterative fixed point control kernel size': np.std(sizesIterativeFilteredFixedPoint),
             'mean iterative rounds': np.mean(iterativeRoundsFiltered),
             'std iterative rounds': np.std(iterativeRoundsFiltered),
             'mean iterative rounds fixed point': np.mean(iterativeRoundsFilteredFixedPoint),
             'std iterative rounds fixed point': np.std(iterativeRoundsFilteredFixedPoint),
             'has limit cycles': d.get('sampled_has_limit_cycles',None),
             'number of attractors': len(d['sampled_attractors']),
             'number of attractors with control kernel': len(attractorsFiltered),
             'number of fixed point attractors': len(attractorsFilteredFixedPoint),
             'all control kernels are identical': ck_identical,
             'iterative control kernel time': d.get('sampled_iterative_control_kernel_time_minutes',None),
             'control kernel time': d.get('sampled_exact_control_kernel_time_minutes',None),
            }
        ddata.update(ddataCK)
    # end control kernel data
    
    return ddata

def loadDataStableMotif(dir='.'):
    dataDict = {}
    dataForFrame = []
    for filename in glob.glob(dir+'/stable_motifs_*.dat'):
        if filename.find('split') == -1: # we don't want to include "split" data
            try:
                d = load(filename)
                success = True
            except:
                print("loadDataStableMotif: Error loading file {}".format(filename))
                success = False
        else:
            success = False
            
        if success:
            # attempting to cure problems with loading python 2 data
            for key in list(d.keys()):
                if type(key) == bytes:
                    d[key.decode()] = d[key]
        
            dataDict[d['name']] = d
            
            ddata = dataFrameStableMotif(d)
            
            if len(ddata) > 0:
                dataForFrame.append(ddata) #dataDictForFrame[d['name']] = ddata

    df = pd.DataFrame.from_records(dataForFrame) #,'index') # from_dict
    df.set_index('name',inplace=True)

    return dataDict,df

def dataFrameStableMotif(d):
    
    ddata = {'name':d['name']}
    
    controlSizes = [ len(ck) for ck in d['control_sets'] ]
    
    ddataSM = {
         'mean control set size': np.mean(controlSizes),
         'std control set size': np.std(controlSizes),
         'number of attractors': len(d['quasi_attractors']),
         }
    ddata.update(ddataSM)
    
    return ddata
    
    

if __name__=='__main__':

    dataDict,df = loadData()
    
    
