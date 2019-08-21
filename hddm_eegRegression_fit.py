
import os
import hddm
from patsy import dmatrix  # for generation of (regression) design matrices
import numpy as np         # for basic matrix operations

model_dir = os.getcwd()

############################################
### 1. FIT THE MODEL #######################
############################################
data = hddm.load_csv("EEGdata_hddm.csv")

samples = 10000 #run 10 models of 10000 trials

#you have to define the specific link functions for z and v 
def z_link_func(x, data=data):
    stim = (np.asarray(dmatrix('0 + C(s,[[1],[-1]])', {'s':data.stimulus.ix[x.index]})))
    return 1 / (1 + np.exp(-(x * stim)))

def v_link_func(x, data=data):
    stim = (np.asarray(dmatrix('0 + C(s,[[1],[-1]])', {'s':data.stimulus.ix[x.index]})))
    return x * stim

#fit different models for previous and following cj
a_reg = {'model': 'a ~ 1 + stimulus + prevpe_bin + postpe_bin + prevernfcz_bin + posternfcz_bin', 'link_func': lambda x: x}
v_reg = {'model': 'v ~ 1 + stimulus + prevpe_bin + postpe_bin + prevernfcz_bin + posternfcz_bin', 'link_func': v_link_func}
reg_descr = [a_reg, v_reg]
m = hddm.HDDMRegressor(data, reg_descr, group_only_regressors=False, p_outlier=.05)
m.find_starting_values()
m.sample(samples, burn=samples/10, thin=2, dbname=os.path.join(model_dir, 'ERPall_traces_1'), db='pickle')
m.save(os.path.join(model_dir, 'ERPall_1'))

goOn = False
if goOn == True:
    import kabuki
    import seaborn as sns
    import matplotlib.pyplot as plt
    models = []
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]: #this assumes you ran 10 different models above and changed the index when saving
        models.append(hddm.load('ERPall_%s' %i))
    m = kabuki.utils.concat_models(models) #contact these 10 models

    #Run the same sanity checks as reported in hddm_fit.py, these are not reiterated here  
    
    #Extract the data and save these
    results = m.gen_stats()
    results.to_csv(os.path.join(model_dir, 'eegRegression_HDDMestimates.csv'))

    vprevern,vprevpe = m.nodes_db.node[['v_prevernfcz_bin', 'v_prevpe_bin']]
    aprevern,aprevpe = m.nodes_db.node[['a_prevernfcz_bin', 'a_prevpe_bin']]
    vpostern,vpostpe = m.nodes_db.node[['v_posternfcz_bin', 'v_postpe_bin']]
    apostern,apostpe = m.nodes_db.node[['a_posternfcz_bin', 'a_postpe_bin']]
    #Save all the traces - and just do the stuff in R!
    numpy.savetxt("v_prevpe_binnedpertrial.csv", vprevpe.trace(), delimiter=",")
    numpy.savetxt("v_prevern_binnedpertrial.csv", vprevern.trace(), delimiter=",")
    numpy.savetxt("a_prevpe_binnedpertrial.csv", aprevpe.trace(), delimiter=",")
    numpy.savetxt("a_prevern_binnedpertrial.csv", aprevern.trace(), delimiter=",")
    numpy.savetxt("v_postpe_binnedpertrial.csv", vpostpe.trace(), delimiter=",")
    numpy.savetxt("v_postern_binnedpertrial.csv", vpostern.trace(), delimiter=",")
    numpy.savetxt("a_postpe_binnedpertrial.csv", apostpe.trace(), delimiter=",")
    numpy.savetxt("a_postern_binnedpertrial.csv", apostern.trace(), delimiter=",")

 
    
