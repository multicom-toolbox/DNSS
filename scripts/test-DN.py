# This script is used to generate and use the DN
# with the specified inputs. It will give the 
# prediction based on the specified model file

# Latest Mod:    11/02/16 by Jie Hou

import sys
GLOBAL_PATH='/space1/jh7x3/DNSS_release/';
sys.path.insert(0, GLOBAL_PATH+'/lib/')
sys.path.insert(0, GLOBAL_PATH+'/lib/cudamat-rbm')

from DN import DN, DN_load,loadmodel 
import numpy as np
import rbm_numpy
import os

if len(sys.argv) < 6:
   sys.stderr.write('Usage: sys.argv[0] test_data_filename model_filename pred probs device <targsize>')
   print "\n"
   sys.exit(1)

DN_model = sys.argv[2]
if not os.path.exists(DN_model):
#if True:
   #dn.test_DeepLearningPro(test_data, target_data)
   print("Error: Model ",DN_model," doesn't exists!\n");
   exit(-1);

device = sys.argv[5]

test_data = np.load(sys.argv[1]) 
test_l1 = test_data[:,:]

if device == 'GPU':
    print "\n!!!Calling GPU for prediction!!!"
    dn = DN_load(DN_model)
    probs = dn.calc_output_legacy(test_l1, 1000)
else:
    print "\n!!!Calling CPU for prediction!!!"
    ###############  load the models 
    loadmodel(DN_model, globals())
    data_l1 = test_l1
    data_l2 = rbm_numpy.calc_hidden_probs(data_l1, l1_vh, l1_hb)
    del data_l1

    data_l3 = rbm_numpy.calc_hidden_probs(data_l2, l2_vh, l2_hb)
    del data_l2

    data_l4 = rbm_numpy.calc_hidden_probs(data_l3, l3_vh, l3_hb)
    del data_l3

    data_l5 = rbm_numpy.calc_hidden_probs(data_l4, l4_vh, l4_hb)
    del data_l4

    data_l6 = rbm_numpy.calc_hidden_probs(data_l5, l5_vh, l5_hb)
    del data_l5
    
    probs = data_l6

targsize = 3
if len(sys.argv) >= 7:
   targsize = int(sys.argv[6])


max_vals = np.reshape(np.repeat(probs.max(axis=1), targsize), (probs.shape[0], targsize))

print "".format(probs[0], max_vals[0], (probs[0] >= max_vals[0]))
preds = 1 * (probs > max_vals - .0001)

np.savetxt(sys.argv[3], preds, fmt="%d")
np.savetxt(sys.argv[4], probs, fmt="%.6f")


