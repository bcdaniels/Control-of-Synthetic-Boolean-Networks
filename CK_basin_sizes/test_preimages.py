import unittest
import preimages as pre
import numpy as np

from neet.boolean.examples import s_pombe, s_cerevisiae
from neet.boolean import LogicNetwork, WTNetwork

# a logic network
from neet.boolean.examples import mouse_cortical_7B
# a weight-threshold network
from neet.boolean.examples import c_elegans

TEST_CELL_COLLECTIVE = True
CELL_COLLECTIVE_DIR = '/Users/bdaniel6/ASUDropbox/Research/control-continuous/Data/Cell Collective/'

def load_cell_collective_network(directory):
    """
    Creates neet network from cell collective data in the
    given directory.
    """
    expressions_file = '{}/expressions.txt'.format(directory)
    external_file = '{}/external.txt'.format(directory)
    net = LogicNetwork.read_logic(expressions_file, external_file)
    return net

def simple_logic_net():
    test_table = [((1,0),{'10','11'}),
                  ((1,),{'1'})]
    return LogicNetwork(test_table)
    
def simple_biased_net(N = 5):
    return WTNetwork(np.ones([N,N]))
    
def simple_dependence_net():
    """
    Nothing depends on node 0
    """
    test_table = [((1,),{'1',}),
                  ((1,),{'1',})]
    return LogicNetwork(test_table)

class TestPreimagesHelpers(unittest.TestCase):
    
    def test_activating_states_sum(self):
        """
        check that the sum of activating and deactivating states
        is 2^(number of inputs)
        """
        for net in [s_pombe,
                    s_cerevisiae,
                    mouse_cortical_7B,
                    c_elegans,
                    simple_logic_net(),
                    simple_biased_net(),
                    simple_dependence_net()]:
            for i in range(net.size):
                neighbors,act_conditions = pre.activating_states(net,i,True)
                neighbors,deact_conditions = pre.activating_states(net,i,False)
                self.assertTrue(
                    len(act_conditions) + len(deact_conditions) == 2**(len(neighbors)))
    
    def test_activating_states_logic(self):
        """
        check that logic table networks have the same inputs
        and activating conditions
        """
        for net in [mouse_cortical_7B,
                    simple_logic_net(),
                    simple_dependence_net()]:
            for i in range(net.size):
                neighbors,act_conditions = pre.activating_states(net,i)
                
                # check that neighbors are the same
                self.assertTrue(tuple(sorted(net.table[i][0])) == neighbors)
                
                # check that length of activating conditions is the same
                self.assertTrue(len(net.table[i][1]) == len(act_conditions))
                
                # check that activating conditions are the same
                if len(neighbors) > 0:
                    table_conditions = [
                        np.array(
                            [int(val) for val in list(condition)]) for condition in net.table[i][1] ]
                    sorted_table_conditions = [
                        tuple(condition[np.argsort(net.table[i][0])]) \
                                    for condition in table_conditions ]
                    self.assertTrue(set(sorted_table_conditions) == set(act_conditions))
            
class TestPreimages(unittest.TestCase):
    
    def test_biased(self):
        for N in [5,10]:
            images = pre.preimages(simple_biased_net(N),np.zeros(N))
            self.assertTrue(len(images),1)
            np.testing.assert_array_equal(np.array(images.iloc[0]),
                             np.zeros(N))
                             
    def test_dimension(self):
        """
        The dimension of preimages should equal the number of nodes.
        """
        for net in [s_pombe,
                    s_cerevisiae,
                    mouse_cortical_7B,
                    c_elegans,
                    simple_logic_net(),
                    simple_biased_net(),
                    simple_dependence_net(),]:
            N = net.size
            images = pre.preimages(net,np.zeros(N))
            self.assertEqual(len(images.columns),N)
    
    def test_update(self):
        """
        Updating the preimages should result in the initial state.
        """
        for net in [s_pombe,
                    s_cerevisiae,
                    mouse_cortical_7B,
                    c_elegans,
                    simple_logic_net(),
                    simple_biased_net(),
                    simple_dependence_net(),]:
            N = net.size
            initial_state = np.zeros(N)
            images = pre.preimages(net,initial_state)
            for i in images.index:
                pre_image_state = list(images.loc[i])
                post_image_state = net.update(pre_image_state)
                np.testing.assert_array_equal(post_image_state,
                                              initial_state)
        
if __name__ == "__main__":
    unittest.main()
