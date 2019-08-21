import os
import hddm
from patsy import dmatrix  # for generation of (regression) design matrices
import numpy as np         # for basic matrix operations

model_dir = os.getcwd()

#########################
### 1. FIT FULL MODEL REGRESSION CODING ###
#########################
### Load Data ###
data = hddm.load_csv("Experiment1_hddm.csv")

samples = 10000 #Run 10 models with 10000 trials and concat the chains

#Define the specific link functions for z and v 
def z_link_func(x, data=data):
    stim = (np.asarray(dmatrix('0 + C(s,[[1],[-1]])', {'s':data.stimulus.ix[x.index]})))
    return 1 / (1 + np.exp(-(x * stim)))

def v_link_func(x, data=data):
    stim = (np.asarray(dmatrix('0 + C(s,[[1],[-1]])', {'s':data.stimulus.ix[x.index]})))
    return x * stim

#### Recode previous and following cj in three bins, correct is the reference
data.prevconf[data.prevconf==1] = 'C_error' #sure error
data.prevconf[data.prevconf==2] = 'C_error' #prob error
data.prevconf[data.prevconf==3] = 'B_guess' #guess error
data.prevconf[data.prevconf==4] = 'B_guess' #guess cor 
data.prevconf[data.prevconf==5] = 'A_correct' #prob cor
data.prevconf[data.prevconf==6] = 'A_correct' #sure cor => this becomes the reference (because of the A)

#### Recode previous and following cj in three bins
data.follconf[data.follconf==1] = 'C_error' #sure error
data.follconf[data.follconf==2] = 'C_error' #prob error
data.follconf[data.follconf==3] = 'B_guess' #guess error
data.follconf[data.follconf==4] = 'B_guess' #guess cor 
data.follconf[data.follconf==5] = 'A_correct' #prob cor
data.follconf[data.follconf==6] = 'A_correct' #sure cor => this becomes the reference

#fit model in which bound an drift depend on prevconf and follconf
a_reg = {'model': 'a ~ 1 + stimulus + prevconf + follconf', 'link_func': lambda x: x}
v_reg = {'model': 'v ~ 1 + stimulus + prevconf + follconf + coh_single', 'link_func': v_link_func}
reg_descr = [a_reg, v_reg]
m = hddm.HDDMRegressor(data, reg_descr, group_only_regressors=True, p_outlier=.05)
m.find_starting_values()
m.sample(samples, burn=samples/10, thin=2, dbname=os.path.join(model_dir, 'Experiment1_traces_1'), db='pickle')
m.save(os.path.join(model_dir, 'Experiment1_1'))

goOn = False
if goOn == True:
    import kabuki
    import seaborn as sns
    import matplotlib.pyplot as plt
    models = []
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]: #this assumes you ran 10 different models above and changed the index when saving
        models.append(hddm.load('Experiment1_%s' %i))
    m = kabuki.utils.concat_models(models) #contact these 10 models
    gelman_rubin(models) #check R hat 

    #diagnostics
    m.plot_posteriors()
    
    #Simulate data and compare to actual data
    data['response'] = data['cor']
    ppc_data = post_pred_gen(m,append_data=True)
    ppc_data['resp_sampled'] = 1
    ppc_data['resp_sampled'][ppc_data.rt_sampled>0] = 0
    hddm.utils.post_pred_stats(data,ppc_data)
        
    ppc_data.rt = abs(ppc_data.rt)
    ppc_data.rt[ppc_data.cor==1] = -ppc_data.rt[ppc_data.cor==1]
    
    sns.distplot(ppc_data.rt_sampled,label='fit',hist=False,kde=True,kde_kws = {'shade': True, 'linewidth': 3}, hist_kws={'edgecolor':'black'})
    sns.distplot(ppc_data.rt,label='data',hist=False,kde=True,kde_kws = {'shade': True, 'linewidth': 3}, hist_kws={'edgecolor':'black'})
    plt.xlim(-5,5)
    
    #Get the results
    results = m.gen_stats()
    results.to_csv(os.path.join(model_dir, 'Experiment1_HDDMestimates.csv'))

    v_guess_prevcj, v_error_prevcj = m.nodes_db.node[['v_prevconf[T.B_guess]', 'v_prevconf[T.C_error]']]
    a_guess_prevcj, a_error_prevcj = m.nodes_db.node[['a_prevconf[T.B_guess]', 'a_prevconf[T.C_error]']]
    v_guess_postcj, v_error_postcj= m.nodes_db.node[['v_follconf[T.B_guess]', 'v_follconf[T.C_error]']]
    a_guess_postcj, a_error_postcj= m.nodes_db.node[['a_follconf[T.B_guess]', 'a_follconf[T.C_error]']]
  
    v_guess = v_guess_prevcj.trace()-v_guess_postcj.trace()
    v_error = v_error_prevcj.trace()-v_error_postcj.trace()
    a_guess = a_guess_prevcj.trace()-a_guess_postcj.trace()
    a_error = a_error_prevcj.trace()-a_error_postcj.trace()
    
    #distributions 
    sns.distplot(a_error,label='perceived error',hist=False,kde=True,kde_kws={'shade': True, 'linewidth': 3},hist_kws={'edgecolor':'black'})
    sns.distplot(a_guess,label='low  confidence',hist=False,kde=True,kde_kws={'shade': True, 'linewidth': 3},hist_kws={'edgecolor':'black'})
    plt.legend();plt.xlabel('Decision bound');plt.ylabel('Posterior distribution')    

    sns.distplot(v_error,label='perceived error',hist=False,kde=True,kde_kws={'shade': True, 'linewidth': 3},hist_kws={'edgecolor':'black'})
    sns.distplot(v_guess,label='low  confidence',hist=False,kde=True,kde_kws={'shade': True, 'linewidth': 3},hist_kws={'edgecolor':'black'})
    plt.legend();plt.xlabel('Drift rate');plt.ylabel('Posterior distribution')    

    #stats
    (a_guess_prevcj.trace()<0).mean()
    (a_error_prevcj.trace()<0).mean()
        
    (a_guess_postcj.trace()>0).mean()
    (a_error_postcj.trace()>0).mean()
    
    (a_guess<0).mean()
    (a_error<0).mean()
    (a_error<a_guess).mean()
    
    (v_guess>0).mean()
    (v_error<0).mean()
    (v_error<v_guess).mean()
    
    #Save these results
    numpy.savetxt("v_guess_prevcj.csv", v_guess_prevcj.trace(), delimiter=",")
    numpy.savetxt("a_guess_prevcj.csv", a_guess_prevcj.trace(), delimiter=",")
    numpy.savetxt("v_error_prevcj.csv", v_error_prevcj.trace(), delimiter=",")
    numpy.savetxt("a_error_prevcj.csv", a_error_prevcj.trace(), delimiter=",")
    numpy.savetxt("v_guess_postcj.csv", v_guess_postcj.trace(), delimiter=",")
    numpy.savetxt("a_guess_postcj.csv", a_guess_postcj.trace(), delimiter=",")
    numpy.savetxt("v_error_postcj.csv", v_error_postcj.trace(), delimiter=",")
    numpy.savetxt("a_error_postcj.csv", a_error_postcj.trace(), delimiter=",")

