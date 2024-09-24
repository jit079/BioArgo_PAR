from pygam.pygam import LinearGAM, s
import numpy as np
import pickle
from pandas import read_csv
import matplotlib.pyplot as plt

def reconstruct_par(profile):
    # profile is a 2-d array, n depth x 5 columns, i.e., [Ed380,Ed443,Ed490,Ed555,depth]
    par=np.zeros((profile.shape[0],100))+np.nan
    for j in range(100):
        model_pkl_file = "gam_3u_2c/bioargo_par_model_"+str(j+1)+".pkl"  
        with open(model_pkl_file, 'rb') as file:  
            model = pickle.load(file)

            par[:,j]=model.predict(profile)

    e1=np.std(par,1)
    
    model_pkl_file = "gam_3u_2c/bioargo_par_uncertainty.pkl"  
    with open(model_pkl_file, 'rb') as file:  
        model = pickle.load(file)  
        e2=model.predict(np.mean(par,1))

    e=np.sqrt(e1**2+e2**2)

    return np.mean(par,1),e

# read example data, 1st column - depth, 2nd - PAR (mW/cm2/micron)
# 3rd - 6th columns: Ed380 (mW/cm2/micron), Ed443 (mW/cm2/micron), Ed490 (mW/cm2/micron), Ed555 (mW/cm2/micron)

df=read_csv('example/example_data.csv',skiprows=None)
df=df.values

d=df[:,0]
par=df[:,1]
profile=df[:,[2,3,4,5,0]]

par_est,uncertainty=reconstruct_par(profile)

print('Estimated PAR:\n')
print(par_est)

plt.figure(figsize=(5,4))
plt.plot(d,par_est,'b',label='modeled')
plt.fill_between(d,par_est-uncertainty,par_est+uncertainty,facecolor='b',alpha=0.2)
plt.plot(d,par,'r',label='measured')
plt.xlabel('Depth, m',fontsize=14)
plt.ylabel('PAR, mW/cm2/micron',fontsize=14)
plt.tick_params(labelsize=12)
plt.yscale('log')
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('example/example_measured_vs_modeled.png')
plt.close()
