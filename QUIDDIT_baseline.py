


###############################################################################
############################ IMPORT SETTINGS #################################

import QUIDDIT_utility as utility
import QUIDDIT_settings as settings
import numpy as np
import os
import scipy.optimize as op
import sys
import matplotlib.pyplot as plt

###############################################################################
############################# INPUT AND DATA ##################################
def main(arg1, arg2):
    #filename=input()
    filename = arg1
    output_path = arg2

#filename='C:\FTIR/aafakedata/LS Arg 81 HQ linescan 40 300001.CSV'
#input_path = 'C:\FTIR\LS Arg 09 HQ linescan'
    #output_path = 'C:\FTIR/aafakedata corrected/'
#if not os.path.exists(IIa_path):
#    print "I can't find the IIa spectrum."
 

    IIa_spec = np.loadtxt(settings.IIa_path, delimiter = ',')

    sumsqu=[]

    i = 1

#if not os.path.exists(input_path):
#    print "input directory doesn't exist"

    #if not os.path.exists(output_path):
    #    print "output directory doesn't exist"
    #    print "creating directory %s" %output_path
    #    os.makedirs(output_path)


#for root, dirs, files in os.walk(input_path):    
#    for name in files:
#        if os.path.splitext(name)[1] == '.CSV' or os.path.splitext(name)[1] == '.csv':
#            spectrumfiles.append(os.path.join(root,name))
#            filenames.append(name)    

#print('reading spectrum %i of %i from: %s...' %(i, len(filenames), name))
    spectrum_prelim = np.loadtxt(filename, delimiter=',')
    spectrum_prelim = utility.spectrum_slice(spectrum_prelim, 675, 4000)
    
    #plt.figure()
    #plt.subplot(2,1,1)
    #plt.plot(spectrum_prelim[:,0], spectrum_prelim[:,1], 'k.', label='original spectrum')
   #plt.legend(loc='best')
    #ax1=plt.gca()
    #ax1.invert_xaxis()
    print('preliminary correction...')
    bl= -spectrum_prelim[-1][1]
   
    spectrum_abs = spectrum_prelim[:,1] + bl              
    spectrum = np.column_stack((spectrum_prelim[:,0], spectrum_abs))
      
    mindiff = (utility.closest(1992.0, spectrum[:,0]))          # return wavenum closest to 1992
    row = np.where(spectrum == mindiff)[0][0]
    factor = 12.3/abs((spectrum[row,1]))                # calculate scaling factor    
    spectrum[:,1] *= factor
    
    #plt.subplot(2,1,2)
    #plt.plot(spectrum[:,0], spectrum[:,1], 'r.', label='corr. spec. prelim.')
          
###############################################################################
################ FITTING AND SUBTRACTING TYPE IIa SPECTRUM #################### 
                                                  
    print('final fit:')                                                    
    two_phonon_left = utility.spectrum_slice(spectrum, 1500,2312)
    two_phonon_right = utility.spectrum_slice(spectrum, 2391, 3000) #3000
    two_phonon_extra = utility.spectrum_slice(spectrum, 3800, 4000)
    two_phonon = np.vstack((two_phonon_left, two_phonon_right, two_phonon_extra))
#two_phonon = spectrum_slice(spectrum, 1500, 2700)
    two_phonon_wav = np.arange(two_phonon[:,0][0], two_phonon[:,0][-1], 0.1)
    two_phonon_ip = utility.inter(spectrum, two_phonon_wav, inttype='linear')          # interpolate slice of spectrum used for fitting    

    IIa_spec_ip = utility.inter(IIa_spec, two_phonon_wav, inttype='linear')            # interpolate relevant area of type IIa spectrum
    IIa_spec_ip_new = utility.inter(IIa_spec, spectrum[:,0:-1], inttype='linear')
    
    IIa_args = (two_phonon_wav, two_phonon_ip, IIa_spec_ip)     # arguments needed for IIa_fit
    IIa_x0 = [(1, 0, 0)]                                    #initial guess of parameters (normf, poly1, poly2)
    IIa_bounds = [(0.0, None),(None, None),(None, None)]         #(min, max)-pairs for parameters 
    IIa_res = op.minimize(utility.IIa, args=IIa_args, x0=IIa_x0, method='L-BFGS-B', bounds=IIa_bounds)
        

    print(IIa_res)
    
    fit_IIa = utility.IIa_fit(IIa_res.x, spectrum[:,0].reshape(len(spectrum[:,0]),1), spectrum[:,1].reshape(len(spectrum[:,1]),1)) 
    sumsqu.append(IIa_res.fun)
    abs_temp = fit_IIa - IIa_spec_ip_new
    
    spec_temp = np.column_stack((spectrum[:,0] , abs_temp))
    
    #f=16
    #if IIa_res.success == False:
    #plt.figure()
    #plt.subplot(2,1,1)
    #plt.title(name)
    #plt.plot(spectrum_prelim[:,0], spectrum_prelim[:,1], 'k.', label='original spectrum')
    #plt.legend(loc='best')
    #ax1=plt.gca()
    #ax1.invert_xaxis()
    #    
    #plt.subplot(2,1,2)
    #plt.plot(spectrum[:,0], spectrum[:,1], 'r.', label='corr. spec. prelim.')
    #plt.plot(spectrum[:,0], fit_IIa, 'k.', label='data fitted to IIa')
    #plt.plot(IIa_spec[:,0], IIa_spec[:,1], 'g-', label='type IIa spectrum')
    #plt.plot(spec_temp[:,0], spec_temp[:,1], 'r.', label='corrected spectrum')
    #plt.plot(spec_temp[:,0],np.polyval(IIa_res.x[1:], spec_temp[:,0]), 'b-', label='baseline')
    #plt.axhline(y=0)
    #plt.legend(loc='best', fontsize=f-8)
        #props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        #plt.text(0.5,0.1,'sum sq: %f' %IIa_res.fun, bbox=props, ha='center', va='bottom', transform = ax1.transAxes, fontsize=12)
    #plt.xlabel(r'$\mathrm{\mathsf{wavenumber\/[cm^{-1}]}}$', fontsize=f)
    #plt.ylabel(r'$\mathrm{\mathsf{absorption\/[cm^{-1}]}}$', fontsize=f)
    #plt.tick_params(axis='both', which='major', labelsize=f-4)
    #ax=plt.gca()
    #ax.invert_xaxis()
    
    #pylab.get_current_fig_manager().window.showMaximized()
    #plt.savefig('FTIR/%s old bl.jpg' %name[-7:-4], dpi=300)
    
    print('saving spectrum after IIa subtraction...')
    #np.savetxt('FTIR/corrected/%s' %'c'+name, spec_temp, delimiter=',') 
#np.savetxt(output_path + 'c' + filename.split('/')[-1], spec_temp, delimiter=',') 
    np.savetxt(output_path + '/c' + filename.split('/')[-1], spec_temp, delimiter=',') 
    
    #output_path = 'FTIR/LS Arg 56 HQ map corrected'
    
    print() 
    print('--------------------------------------------------------------------')
    print()

###############################################################################
################################ PLOTTING #####################################    

    #plt.figure()
    #plt.subplot(2,1,1)
    #plt.plot(spectrum_prelim[:,0], spectrum_prelim[:,1], 'k.', label=' original data')
    #plt.tick_params(axis='both', which='major', labelsize=f)
    #plt.xlim(675,4000)
    #ax1=plt.gca()
    #ax1.invert_xaxis()
    #ax1.set_xticks([])
    #plt.ylabel('absorbance ($\mathregular{cm^{-1}}$)', fontsize=f+2)

    #plt.subplot(2,1,2)        
    #plt.plot(spectrum[:,0], spectrum[:,1], 'k.', label='data')
    #plt.plot(spectrum[:,0], fit_IIa, '.', color='C1', label='data fitted to IIa')
    
    #plt.plot(IIa_spec[:,0], IIa_spec[:,1], 'k-', label='standard')
    
    #plt.plot(spec_temp[:,0], spec_temp[:,1], 'g.', label='corrected spectrum')
    #plt.axhline(y=0, linestyle='--', color='0.5', lw=2)
    #plt.plot(spectrum[:,0], np.polyval(IIa_res.x[1:], spectrum[:,0]), '--', color='0.5', label='baseline')

    #plt.legend(loc='upper left', fontsize=f)
    #plt.title(name)
    #plt.suptitle(IIa_res.success)
    #plt.tick_params(axis='both', which='major', labelsize=f)
    #plt.xlim(675,4000)
    #ax2=plt.gca()
    #ax2.invert_xaxis()
    #plt.ylabel('absorbance ($\mathregular{cm^{-1}}$)', fontsize=f+2)
    #plt.xlabel('wavenumber ($\mathregular{cm^{-1}}$)', fontsize=f+2)
    i +=1

   
    #plt.show()

    print('************************************************************************')
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])