###############################################################################
######################### IMPORTS AND FUNCTIONS ###############################
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
from scipy import interpolate
import sys

def inter(spec, wav_new, inttype='cubic'):
    """returns interpolated spectrum"""
    spec_interp = interpolate.interp1d(spec[:,0],spec[:,1], kind=inttype, bounds_error=False, fill_value=0)
    return spec_interp(wav_new) 
				
def closest(target, collection):
    """returns value from collection that is closest to target"""
    return min((abs(target-i),i) for i in collection)[1]

def aggregate(T1, NT, T2, IaB_measured, t1, t2):
    """input: NT,t1,t2, measured agg (prop. of B), T2.
    vary T1, until IaB2 matches IaBt
    """
    EaR = 81160
    preexp = 293608
    NA0 = NT    #A centres before aggregation starts
    rate1 = preexp * np.exp(-EaR/(T1+273))
    NA1 =  NA0/(1+(rate1*NA0)*t1)     #A centres after first stage of annealing
    rate2 = preexp * np.exp(-EaR/(T2+273))
    NA2 = NA1/(1+(rate2*NA1)*t2)      #A centres after second stage of annealing
    IaB_calc = 1-(NA2/NT)
    error = IaB_measured - IaB_calc
    return error**2
    
def temperature(t, NT, IaB):
    NA = NT * (1-IaB)
    T = (-81160/(np.log(((NT/NA)-1)/(t*NT*293608))))
    return T-273

###############################################################################
######################### CHANGES THESE VALUES ################################
def main(arg1, arg2, arg3, arg4, arg5):

    output = 'output.CSV'
# total duration in Ma:
# (both stages together!)
    #age_tot = 530  
    age_tot = arg1
# Nitrogen data:
# core:
    #c_NT =  643         #total N
    #c_agg =  0.515      #proportion of B
    c_NT = arg2
    c_agg = arg3

#rim:
    #r_NT = 85
    #r_agg = 0.045
    r_NT = arg4
    r_agg = arg5

# settings for diagram:
    #f=16                    #fontsize
    #T_limit = (1100, 1500)  #(min, max) values used for T

###############################################################################
########################### SETTING UP THE MODEL ##############################

    stage1 = np.arange(0, age_tot+10, 10)   #possible durations of stage 1 in 10 Ma steps
    stage1[0]=1
    stage1[-1]=age_tot-1
    #T_bound = [(900,1400)]

    #plt.figure()

    st1 = []
    st2 = []


    for duration in stage1:
        r_duration = (age_tot-duration) * 1e6 * 365.25 * 24 * 60 * 60
        c_duration = duration * 1e6 * 365.25 * 24 * 60 * 60
        r_T = temperature(r_duration, r_NT, r_agg)
        st2.append(r_T)
        print("rim T: {}".format(r_T))
        #x0 = 1165
        res = op.minimize_scalar(aggregate, bounds=(700,1500), args=(c_NT,r_T,c_agg,c_duration,r_duration), method='bounded')

        print(res)
        c_T = res.x
        st1.append(c_T)
        #plt.plot(duration, c_T,  'ro')
        #plt.plot(duration, r_T,  'bo')
    
#st1_inter = inter(np.column_stack((stage1,st1)), np.arange(0, age_tot, 2))
#st2_inter = inter(np.column_stack((stage1,st2)), np.arange(0, age_tot, 2))

    #plt.gca().tick_params(labelsize=f-2)
			
    #plt.xlabel('Duration of first anneal (Ma)', fontsize=f)    
    #plt.ylabel('Temperature ($\mathregular{^{\circ}}$C)', fontsize=f)
    #plt.xlim(0, age_tot)
    #plt.ylim(T_limit)
    #plt.tight_layout()
    #plt.show()
    
    np.savetxt(output, np.column_stack((stage1, st1, st2)), delimiter=',', header='duration of first anneal, T_core, T_rim')
    #with open(output, 'w') as output_fob:
        #output_fob.write('duration of first anneal, T_core, T_rim, \n')
    
        #for row in zip(stage1, st1, st2):
            #for item in row:
                #output_fob.write(str(item)+',')
            #output_fob.write('\n')

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
