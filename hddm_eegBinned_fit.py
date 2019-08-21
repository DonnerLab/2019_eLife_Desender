import os
import hddm
from patsy import dmatrix  # for generation of (regression) design matrices
import numpy as np         # for basic matrix operations

model_dir = os.getcwd()

############################################
### 1. FIT THE MODEL #######################
############################################
data = hddm.load_csv("EEGdata_hddm.csv")

samples = 10000 #run 10 models of 1000O

data.prevpe_bin_noern[data.prevpe_bin_noern==1] = 'A_bin' #this becomes the reference
data.prevpe_bin_noern[data.prevpe_bin_noern==2] = 'B_bin' #
data.prevpe_bin_noern[data.prevpe_bin_noern==3] = 'C_bin' #
data.prevpe_bin_noern[data.prevpe_bin_noern==4] = 'D_bin' #
data.prevpe_bin_noern[data.prevpe_bin_noern==5] = 'E_bin' #

data.postpe_bin_noern[data.postpe_bin_noern==1] = 'A_bin' #this becomes the reference
data.postpe_bin_noern[data.postpe_bin_noern==2] = 'B_bin' #
data.postpe_bin_noern[data.postpe_bin_noern==3] = 'C_bin' #
data.postpe_bin_noern[data.postpe_bin_noern==4] = 'D_bin' #
data.postpe_bin_noern[data.postpe_bin_noern==5] = 'E_bin' #

#define the specific link functions for z and v 
def z_link_func(x, data=data):
    stim = (np.asarray(dmatrix('0 + C(s,[[1],[-1]])', {'s':data.stimulus.ix[x.index]})))
    return 1 / (1 + np.exp(-(x * stim)))

def v_link_func(x, data=data):
    stim = (np.asarray(dmatrix('0 + C(s,[[1],[-1]])', {'s':data.stimulus.ix[x.index]})))
    return x * stim

#fit different models for previous and following cj
a_reg = {'model': 'a ~ 1 + stimulus + prevpe_bin_noern + postpe_bin_noern', 'link_func': lambda x: x}
v_reg = {'model': 'v ~ 1 + stimulus + prevpe_bin_noern + postpe_bin_noern', 'link_func': v_link_func}
reg_descr = [a_reg, v_reg]
m = hddm.HDDMRegressor(data, reg_descr, group_only_regressors=False, p_outlier=.05)
m.find_starting_values()
m.sample(samples, burn=samples/10, thin=2, dbname=os.path.join(model_dir, 'ERPall_binnedPE_traces_1'), db='pickle')
m.save(os.path.join(model_dir, 'ERPall_binnedPE_1'))

goOn = False
if goOn == True:
    import kabuki
    import seaborn as sns
    import matplotlib.pyplot as plt
    models = []
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]: #this assumes you ran 10 different models above and changed the index when saving
        models.append(hddm.load('ERPall_binnedPE_%s' %i))
    m = kabuki.utils.concat_models(models) #contact these 10 models

    #Run the same sanity checks as reported in hddm_fit.py, these are not reiterated here  
    
    #Extract the data and save these
    results = m.gen_stats()
    results.to_csv(os.path.join(model_dir, 'eeg_binnedPe_HDDMestimates.csv'))
    
    a2, a3, a4, a5 = m.nodes_db.node[['a_prevpe_bin_noern[T.B_bin]', 'a_prevpe_bin_noern[T.C_bin]', 'a_prevpe_bin_noern[T.D_bin]', 'a_prevpe_bin_noern[T.E_bin]']]
    p2, p3, p4, p5 = m.nodes_db.node[['a_postpe_bin_noern[T.B_bin]', 'a_postpe_bin_noern[T.C_bin]', 'a_postpe_bin_noern[T.D_bin]', 'a_postpe_bin_noern[T.E_bin]']]
    v2, v3, v4, v5 = m.nodes_db.node[['v_prevpe_bin_noern[T.B_bin]', 'v_prevpe_bin_noern[T.C_bin]', 'v_prevpe_bin_noern[T.D_bin]', 'v_prevpe_bin_noern[T.E_bin]']]
    vp2, vp3, vp4, vp5 = m.nodes_db.node[['v_postpe_bin_noern[T.B_bin]', 'v_postpe_bin_noern[T.C_bin]', 'v_postpe_bin_noern[T.D_bin]', 'v_postpe_bin_noern[T.E_bin]']]

    numpy.savetxt("PEonly_noern_a_prevcj_bin2.csv", a2.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_a_prevcj_bin3.csv", a3.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_a_prevcj_bin4.csv", a4.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_a_prevcj_bin5.csv", a5.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_a_postcj_bin2.csv", p2.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_a_postcj_bin3.csv", p3.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_a_postcj_bin4.csv", p4.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_a_postcj_bin5.csv", p5.trace(), delimiter=",")

    numpy.savetxt("PEonly_noern_v_prevcj_bin2.csv", v2.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_v_prevcj_bin3.csv", v3.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_v_prevcj_bin4.csv", v4.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_v_prevcj_bin5.csv", v5.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_v_postcj_bin2.csv", vp2.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_v_postcj_bin3.csv", vp3.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_v_postcj_bin4.csv", vp4.trace(), delimiter=",")
    numpy.savetxt("PEonly_noern_v_postcj_bin5.csv", vp5.trace(), delimiter=",")

