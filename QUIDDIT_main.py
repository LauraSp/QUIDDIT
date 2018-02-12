#-------------------------------------------------------------------------------
# Name:         QUIDDIT main
# Purpose:      spectral deconvolution
#
# Author:      Laura Speich
#
# Created:     20/11/2014
# Copyright:   (c) ls13943 2014
#-------------------------------------------------------------------------------

###############################################################################
####################### IMPORTING REQUIRED MODULES ############################
 
import QUIDDIT_utility as utility
import QUIDDIT_settings as settings
import numpy as np
import sys
import scipy.optimize as op
from scipy import integrate
from scipy import stats
       
###############################################################################
########################### INPUT & DATA ######################################

def main(arg1, arg2, arg3):
    filename = arg1
    age = arg2
    
    #N_comp = arg3
    
    N_selection = np.array((arg3))
    #print(N_selection)

    #N_selection = arg3    


    results = np.zeros(1, dtype=utility.results_dtype)
    review = np.zeros(1, dtype=utility.review_dtype)


    C = np.column_stack((settings.std[:,0], settings.std[:,1]))
    A = np.column_stack((settings.std[:,0], settings.std[:,2]))    #generate C, A, X, B and D std
    X = np.column_stack((settings.std[:,0], settings.std[:,3]))    #spectra from CAXBD file
    B = np.column_stack((settings.std[:,0], settings.std[:,4]))   
    D = np.column_stack((settings.std[:,0], settings.std[:,5]))

###############################################################################
########################### READ AND PROCESS DATA #############################

    spectrum = np.loadtxt(filename, delimiter=',')          # generate np array from file
  
  
###############################################################################
############################ 3107cm-1 HYDROGEN PEAK ###########################

    print('fitting Pseudovoigt function to 3107 cm-1 hydrogen peak...')
    
# extract area around 3107 cm-1 H peak, fit and subtract baseline: 
    H_area = utility.spectrum_slice(spectrum, 3000, 3200)
    H_bg_left = utility.spectrum_slice(spectrum, 3000, 3050)
    H_bg_right = utility.spectrum_slice(spectrum, 3150, 3200)
    H_bg_both = np.vstack((H_bg_left, H_bg_right))
    
    H_p_bg = np.polyfit(H_bg_both[:,0], H_bg_both[:,1], 3)      # fit 3rd order polynomial baseline
    H_bg = np.polyval(H_p_bg, H_area[:,0])

    H_absorp_new = H_area[:,1] - H_bg
    H_wav_new = H_area[:,0]
    spec_new = np.column_stack((H_wav_new, H_absorp_new))
    
    H_bg_a, H_bg_b, H_bg_c, H_bg_d = H_p_bg
    #H_bg_a = H_p_bg[0]
    #H_bg_b = H_p_bg[1] 
    #H_bg_c = H_p_bg[2]
    #H_bg_d = H_p_bg[3]    

# fit Pseudovoigt function to 3107 cm-1 H peak:       
    H_x0 = (3107, 0, 1, 1, 0.5) #first guess for x0, I, HWHM_left and HWHM_right, sigma
    H_bounds = [(3106,3108),(0,None),(0.001,5),(0.001,5),(0,1)]

    wav_inter = np.arange(H_wav_new[0], H_wav_new[-1], 0.1)
    peak_inter = utility.inter(spec_new, wav_inter)

    H_res = op.minimize(utility.pseudovoigt, method='SLSQP', args=(wav_inter, peak_inter), x0=H_x0, bounds=H_bounds)

    H_fit = utility.pseudovoigt_fit(H_wav_new, *H_res.x)
    
    H_pos, H_I, H_HWHM_l, H_HWHM_r, H_sigma = H_res.x
    #H_pos = H_res.x[0]
    #H_I = H_res.x[1]
    #H_HWHM_l.append(H_res.x[2])
    #H_HWHM_r.append(H_res.x[3])
    #H_sigma.append(H_res.x[4])
    
# calculate peak area:
    print('calculating peak area...')          
    H_spec = utility.spectrum_slice(spectrum, 3102, 3112)
    H_spec2 = H_spec[:,1] - np.polyval(H_p_bg, H_spec[:,0])
    #H_area_numerical_data = integrate.simps(H_spec2)    
    #H_area_numerical_fit = integrate.simps(H_fit)      #integrate bg corrected fit
    H_area_analytical = 2*(H_res.x[1])*((H_res.x[2]+H_res.x[3])/2)*(H_res.x[4]*(np.pi/2)+(1-H_res.x[4])*np.sqrt(np.pi/2))      
          
###############################################################################
############################## PLATELET PEAK ##################################   
    
# extract area around platelet peak, fit and subtract preliminary background   
    print('fitting pseudovoigt function to platelet peak...')  
    pp = utility.spectrum_slice(spectrum, 1327, 1420)         # extract area around pp
    pp2 = utility.spectrum_slice(spectrum, 1350, 1380)
    
    #p_s2n = stats.signaltonoise(pp[:,1])
    
    pp_wav_inter = np.arange(pp[0][0], pp[-1][0], 0.01)
    pp_inter = utility.inter(pp, pp_wav_inter)

# calculations for bounds
    I1405 = 0.257 * H_res.x[1]      #from empirical analysis on platelet degraded spectra 
    H_lb = I1405 - .1*I1405
    H_ub = I1405 + .1*I1405

    I1332 = utility.height(1332, pp)
    if I1332 <= 0:
         B_lb = 0
         B_ub = .5
    else:
        B_lb = I1332 - .1*I1332
        B_ub = I1332 + .1*I1332
    
    p_max = pp2[np.argmax(pp2[:,1])][0]
    p_lb = p_max - 1.5
    p_ub = p_max + 1.5
    
    cc = min(pp[:,1])
    if cc>=0:
        c_ub = cc
        c_lb = 0
    else:
        c_ub = 0
        c_lb = cc

# fit pseudovogt functions to platelet peak, 1405 and 1332 at the same time
    psv_x0 = (p_max, 0, 5, 5, 1, 1405, I1405, 5, 5, 1, 1332, I1332, 5, 5, 0, 0)
    p_bounds = [(p_lb,p_ub),(0,None),(.01,50),(.01,50),(0,1),(1404.5,1405.5),(H_lb,H_ub),(.1,5),(.1,5),(0,1),(1331,1333),(B_lb,B_ub),(.1,5),(.1,5),(0,1), (c_lb,c_ub)] 
    cons = ({'type': 'ineq', 'fun': utility.pp_cons1}, {'type': 'ineq', 'fun': utility.pp_cons2}, {'type': 'ineq', 'fun': utility.pp_cons3})    

    pp_res = op.minimize(utility.ultimatepsv, x0=psv_x0, args=(pp_wav_inter, pp_inter), method='SLSQP', bounds=p_bounds, constraints=cons)
    #pp_fit = utility.ultimatepsv_fit(pp_wav_inter, *pp_res.x)                        
    
    if pp_res.success != True:
        print('trying alternative method for fitting platelet peak')
        psv_x0 = settings.pp_res_prev
        pp_res = op.minimize(utility.ultimatepsv, x0=psv_x0, args=(pp_wav_inter, pp_inter), method='SLSQP', bounds=p_bounds, constraints=cons)
        #pp_fit = utility.ultimatepsv_fit(pp_wav_inter, *pp_res.x)                              
        
        
    pp_sumsqu = pp_res.fun     
      
    
    #if pp_res.success == True and pp_res.x[1]>=1:
    if pp_res.x[1]>0:
        p_x0, p_I, p_HWHM_l, p_HWHM_r, p_sigma, B_x0, B_I, B_HWHM_l, B_HWHM_r, B_sigma, H1405_x0, H1405_I, H1405_HWHM_l, H1405_HWHM_r, H1405_sigma, psv_c = pp_res.x
        #p_x0.append(pp_res.x[0])
        #p_I.append(pp_res.x[1])
        #p_HWHM_l.append(pp_res.x[2])
        #p_HWHM_r.append(pp_res.x[3])
        #p_sigma.append(pp_res.x[4])
        #B_x0.append(pp_res.x[5])
        #B_I.append(pp_res.x[6])
        #B_HWHM_l.append(pp_res.x[7])
        #B_HWHM_r.append(pp_res.x[8])
        #B_sigma.append(pp_res.x[9])
        #H1405_x0.append(pp_res.x[10])
        #H1405_I.append(pp_res.x[11])
        #H1405_HWHM_l.append(pp_res.x[12])
        #H1405_HWHM_r.append(pp_res.x[13])
        #H1405_sigma.append(pp_res.x[14])
        #psv_c.append(pp_res.x[15])
        
        
        # calculate peak area in different ways:
        print('calculating peak area and peak symmetry...')
        if pp_res.x[2]<1:
            int_lower_bound = pp_res.x[0] - 15
        else:
            int_lower_bound = pp_res.x[0] - 15*pp_res.x[2]
        
        if pp_res.x[3]<1:
            int_upper_bound = pp_res.x[0] + 15
        else:
            int_upper_bound = pp_res.x[0] + 15*pp_res.x[3]
            
        pp_spec = utility.spectrum_slice(spectrum, int_lower_bound, int_upper_bound)
        
        #H1405_psv = utility.pseudovoigt_fit(pp_spec[:,0],*pp_res.x[5:10])
        
        #B_psv = utility.pseudovoigt_fit(pp_spec[:,0],*pp_res.x[10:-1])
        
        pp_abs_new = utility.pseudovoigt_fit(pp_spec[:,0],*pp_res.x[:5])
            
        p_area_numerical_data = integrate.simps(pp_abs_new)
        #p_area_numerical_fit = integrate.simps(utility.pseudovoigt_fit(pp_spec[:,0],*pp_res.x[:5]))      #integrate bg corrected fit
        p_area_analytical = 2*(pp_res.x[1])*((pp_res.x[2]+pp_res.x[3])/2)*(pp_res.x[4]*(np.pi/2)+(1-pp_res.x[4])*np.sqrt(np.pi/2))                                       
                                                                                                                                        

# symmetry calculations:
        sym_lower_bound = pp_res.x[0]-5*pp_res.x[2]
        sym_upper_bound = pp_res.x[0]+5*pp_res.x[3]
        
        wavenum_l = np.linspace(sym_lower_bound, pp_res.x[0], 100)
        wavenum_r = np.linspace(pp_res.x[0], sym_upper_bound, 100)
        
        absorp_l = utility.pseudovoigt_fit(wavenum_l, *pp_res.x[:5])
        absorp_r = utility.pseudovoigt_fit(wavenum_r, *pp_res.x[:5])         

# peak asymmetry factor (As):    
        I10 = 0.1 * pp_res.x[1]    
        wav10_l=np.where(absorp_l == utility.closest(I10, absorp_l))[0][0]
        wav10_r=np.where(absorp_r == utility.closest(I10, absorp_r))[0][0]

        a_As = wavenum_r[wav10_r]-pp_res.x[0]
        b_As = pp_res.x[0] - wavenum_l[wav10_l]
    
        asymmetry_factor=b_As/a_As

    
# tailing factor (Tf):
        I5 = 0.05 * pp_res.x[1]
        wav5_l=np.where(absorp_l == utility.closest(I5, absorp_l))[0][0]
        wav5_r=np.where(absorp_r == utility.closest(I5, absorp_r))[0][0]   

        a_Tf = wavenum_r[wav5_r]-pp_res.x[0]
        b_Tf = pp_res.x[0] - wavenum_l[wav5_l]
            
        tailing_factor=a_Tf+b_Tf/(2*a_Tf)

# integral breadth (beta):
        integral_breadth = p_area_numerical_data/pp_res.x[1]
    
# form factor (phi):
        form_factor = (pp_res.x[2]+pp_res.x[3])/integral_breadth
    
# centroid/weighted average:
        avg=(np.average(pp_spec[:,0], weights=pp_abs_new))         
        
        
    else:
        print('no platelet peak found')
        
        
        p_x0 = p_I = p_HWHM_l = p_HWHM_r = p_sigma = B_x0 = B_I = B_HWHM_l = B_HWHM_r = B_sigma = H1405_x0 = H1405_I = H1405_HWHM_l = H1405_HWHM_r = H1405_sigma = psv_c = np.nan
        #p_x0.append(np.nan)
        #p_I.append(np.nan)
        #p_HWHM_l.append(np.nan)
        #p_HWHM_r.append(np.nan)
        #p_sigma.append(np.nan)
        #B_x0.append(np.nan)
        #B_I.append(np.nan)
        #B_HWHM_l.append(np.nan)
        #B_HWHM_r.append(np.nan)
        #B_sigma.append(np.nan)
        #H1405_x0.append(np.nan)
        #H1405_I.append(np.nan)
        #H1405_HWHM_l.append(np.nan)
        #H1405_HWHM_r.append(np.nan)
        #H1405_sigma.append(np.nan) 
        #psv_c.append(np.nan)       
        
        p_area_numerical_data = p_area_numerical_fit = p_area_analytical = asymmetry_factor = tailing_factor = integral_breadth = form_factor = avg = np.nan                                           

    pp_res_prev = pp_res.x

###############################################################################
########################### NITROGEN AGGREGATION ##############################        
    
# interpolate A, B and D spectra and sample spectrum:
    print('calculating N aggregation data...')
    N_area = utility.spectrum_slice(spectrum, 1001, 1399)
    N_wav_new = N_area[:,:-1]  

    C_new = utility.inter(C, N_wav_new)   
    A_new = utility.inter(A, N_wav_new)
    X_new = utility.inter(X, N_wav_new)                     # interpolate C, A, X, B and D spectra
    B_new = utility.inter(B, N_wav_new)
    D_new = utility.inter(D, N_wav_new)
    
    N_area_new = utility.inter(N_area, N_wav_new)          # interpolate measured spectrum (N part)
    N_area_new = N_area_new.flatten()
    
# initial parameters and bounds for fit:
    min_area = 5
    
    if p_area_analytical < min_area:
        N_x0_d = 0.0
    else:
        N_x0_d = None        
    
    #polyx0 = N_area[-1,1]    
    #if polyx0 >0:        
    #    polybounds = (0, polyx0)    
    #else:
    #    polybounds = (polyx0, 0)
    
    if N_selection[-1] == 1:    
        
        polyx0 = N_area[-1,1]    
        if polyx0 >0:        
            polybounds = (0., polyx0)    
        else:
            polybounds = (polyx0, 0.)
    else:
        polyx0 = 0
        polybounds = (0.,0.)
    
                
    #N_x0 = (.5, .5, .1, .5, 0., -polyx0)[np.where(N_selection==1)]                     # initial guess for a, b, d and poly1
    N_x0 = [i for i,j in zip((.5, .5, .1, .5, 0., -polyx0), N_selection) if j==1]
    #N_x0=[(.5, .5, .5, 0, -polyx0)]                               # initial guess for c, a, b, d and poly1
    #N_bounds = [(0.,None),(0.,None),(0.,None),(0.,None),(0., N_x0_d), polybounds][np.where(N_selection==1)]
    N_bounds =  [i for i,j in zip([(0.,None),(0.,None),(0.,None),(0.,None),(0., N_x0_d), polybounds], N_selection) if j==1]
    #N_bounds = [(0.0, None),(0.0,None),(0.0,None),(0.0, N_x0_d), polybounds]   # bounds for c, a, b and d ((min, max)-pairs)
    #N_bounds = [(0.0,None),(0.0,None),(0.0, N_x0_d), (0,0)]   # bounds for a, b and d ((min, max)-pairs)

    N_args = np.column_stack((C_new, A_new, X_new, B_new, D_new))[:,np.where(N_selection[:-1]==1)[0]]     # arguments needed for least-squares fit of ABD-function
    #N_args = (N_area_new, C_new, A_new, B_new, D_new)
    
    if N_selection[-2]==1:
        if N_selection[0] == 1:
            N_cons = ({'type': 'ineq', 'fun': utility.Nd_bound2})
        else:
            N_cons = ({'type': 'ineq', 'fun': utility.Nd_bound})
        
    else:
        N_cons = None



# optimization: 
    N_res = op.minimize(utility.CAXBD_err, x0=N_x0, args=(N_args, N_area_new), method='SLSQP', bounds=N_bounds, constraints=N_cons)
    #N_res = op.minimize(CABD, x0=N_x0, args=N_args, method='SLSQP', bounds=N_bounds, constraints=N_cons)
    
    print(N_res)
    
    c, a, x, b, d, N_poly = [], [], [], [], [], []
    #c, a, x, b, d, N_poly = None, None, None, None, None, None
    N_res_idx = 0
    
    for (idx, val), comp in zip(enumerate(N_selection), (c, a, x, b, d, N_poly)):
        if val == 0:
            comp.append(np.nan)
        elif val == 1:
            comp.append(N_res.x[N_res_idx])
            N_res_idx += 1
            

    N_sumsqu = N_res.fun     # sum of squares of measured - fit

# fit to data (for plotting):
    #N_fit = utility.ABD_fit(N_res.x[0], N_res.x[1], N_res.x[2], N_res.x[3], A_new, B_new, D_new)
    #N_fit = CABD_fit(N_res.x[0], N_res.x[1], N_res.x[2], N_res.x[3], N_res.x[4], C_new, A_new, B_new, D_new)
    
    #plt.figure()
    #plt.plot(N_wav_new, N_fit, label='ABD fit')
    ##plt.plot(N_wav_new, N_res.x[0]*C_new, label='C')
    #plt.plot(N_wav_new, N_res.x[1]*A_new, label='A')
    #plt.plot(N_wav_new, N_res.x[2]*B_new, label='B')
    #plt.plot(N_wav_new, N_res.x[3]*D_new, label='D')
    ##plt.axhline(y=N_res.x[3])
    #plt.plot(N_area[:,0], N_area[:,1], '.', label='data')
    #plt.legend(loc='best')
    #ax=plt.gca()
    #ax.invert_xaxis()

# calculate concentrations of A- and B-centres and total N (in ppm):    
    NC = float(c[0])*25
    NA = float(a[0])*16.5
    NB = float(b[0])*79.4
    NT = np.sum(np.nan_to_num((NC, NA, NB)))
    
    #perc_IaB = NB/(NA+NB)        # %aggregation
    #Woods = 50*float(b[0])    # expected platelet peak area (Woods 1986)
    
###############################################################################
############################### TEMPERATURE ###################################    

# note: in python "ln" is log
    age_s = age * 1e6 * 365 * 24 * 60 * 60                  #age in seconds
    T = (-81160/(np.log(((NT/NA)-1)/(age_s*NT*293608)))) - 273.15    
    
    if T < 0:
        temperature= np.nan
    else:
        temperature = T         
           
                     
###############################################################################
##################### WRITE TO RESULTS FILE ###################################
																																									
    print('saving results to file...')
    #res_header = 'name, p_x0, p_I, p_HWHM_l, p_HWHM_r, p_sigma, avg, area_num_data, area_num_fit, area_ana, As, Tf, beta, phi, pp_sumsqu, p_s2n, c, a, b, d, const, [NA], [NB], [Nt], perc_IaB, Woods, T, N_sumsqu, I_3107, H_area_num_data, H_area_num_fit, H_area_ana'
          
# write all results to structured arrays and store them as csv-files:                          
    
    results['name'] = filename

    results['p_x0'] = p_x0
    results['p_I'] = p_I
    results['p_HWHM_l'] = p_HWHM_l
    results['p_HWHM_r'] = p_HWHM_r
    results['p_sigma'] = p_sigma
    results['p_area_num_data'] = p_area_numerical_data
    results['p_area_ana'] = p_area_analytical
    results['avg'] = avg
    results['p_As'] = asymmetry_factor
    results['p_Tf'] = tailing_factor
    results['p_beta'] = integral_breadth
    results['p_phi']  = form_factor
    results['p_sumsqu'] = pp_sumsqu

    results['c'] = c
    results['a'] = a  
    results['x'] = x       
    results['b'] = b         
    results['d'] = d
    results['N_poly'] = N_poly
    results['[NC]'] = NC
    results['[NA]'] = NA        
    results['[NB]'] = NB                 
    results['[NT]'] = NT
    results['T'] = temperature
    results['N_sumsqu'] = N_sumsqu

    results['I_3107'] = H_I
    results['H_area_ana'] = H_area_analytical 
												
   
    review['name'] = filename

    review['p_x0'] = p_x0
    review['p_I'] = p_I
    review['p_HWHM_l'] = p_HWHM_l
    review['p_HWHM_r'] = p_HWHM_r
    review['p_sigma'] = p_sigma
    review['B_x0'] = B_x0
    review['B_I'] = B_I
    review['B_HWHM_l'] = B_HWHM_l
    review['B_HWHM_r'] = B_HWHM_r
    review['B_sigma'] = B_sigma
    review['H1405_x0'] = H1405_x0
    review['H1405_I'] = H1405_I
    review['H1405_HWHM_l'] = H1405_HWHM_l
    review['H1405_HWHM_r'] = H1405_HWHM_r
    review['H1405_sigma'] = H1405_sigma
    review['psv_c'] = psv_c
    review['avg'] = avg

    review['c'] = c
    review['a'] = a
    review['x'] = x
    review['b'] = b
    review['d'] = d
    review['N_poly'] = N_poly

    review['H_bg_a'] = H_bg_a
    review['H_bg_b'] = H_bg_b
    review['H_bg_c'] = H_bg_c
    review['H_bg_d'] = H_bg_d
    review['H_pos'] = H_pos
    review['H_I'] = H_I
    review['H_HWHM_l'] = H_HWHM_l
    review['H_HWHM_r'] = H_HWHM_r
    review['H_sigma'] = H_sigma

    review['path'] = filename                         

    return results, review             
                    

###############################################################################
############################## OVERALL PLOTTING ###############################

#ax=plt.gca()
#ax.invert_xaxis()
#plt.title(sample)
#plt.legend(loc='best')
#plt.xlabel(r'$\mathrm{\mathsf{wavenumber\/[cm^{-1}]}}$')
#plt.ylabel(r'$\mathrm{\mathsf{absorption\/[cm^{-1}]}}$')
#plt.show()

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1], sys.argv[2])