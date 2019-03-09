#-*- coding: utf-8 -*-
"""
Created on Thu Sep 07 12:31:42 2017

@author: ls13943
"""

import sys
import os
import webbrowser
import tkinter as tk
import tkinter.filedialog as fd
from tkinter import simpledialog as sd
from tkinter import ttk
import numpy as np
from scipy import interpolate
import scipy.optimize as op

#import spectral as sp
#import spectral.io.envi as envi

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.cm as cm
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider
from matplotlib.widgets import Button as mplButton

from tclwinbase_QUIDDITversion import *
import QUIDDIT_settings as settings
import QUIDDIT_utility as utility
import QUIDDIT_baseline
import QUIDDIT_main
import QUIDDIT_two_stage_model
import QUIDDIT_peakfit

QUIDDITversion = settings.version
STDBG = '#ececec'
STDCOLS = cm.jet

#all_toolitems = ('Home', 'Back', 'Forward', None,
#                 'Pan', 'Zoom', 'Subplots', None,
#                 'Save')


###############################################################################

class QUIDDITMain(TclWinBase):
    """Main window for QUIDDIT
    """

    def make_gui(self, title):
        """
        """
        self.setwintitle(title)

        self.set_defaults()

        self.master.protocol('WM_DELETE_WINDOW', self.click_exit)
        self.master.bind('<Return>', self.click_ok)
        self.master.bind('<Escape>', self.click_exit)

        #Making the main menu
        menubar = tk.Menu(self)

        filemenuopt = {'Open...':self.file_open,
                       'Convert ENVI file...':self.convert_ENVI,
                       'Exit':self.click_exit}
        filemenu = self.make_menu(menubar, 'File', filemenuopt)
        filemenu.insert_separator(2)

        optmenuopt = {'settings':self.settings,
                      'Custom baseline':self.hello}
        self.make_menu(menubar, 'Settings', optmenuopt)

        #baselinemenuopt = {'Correct baseline (default)':self.baseline,
        #                   'Correct with custom baseline':self.hello,
        #                   'Another option':self.hello}
        baselinemenuopt = {'Correct baseline':self.baseline}
        self.make_menu(menubar, 'Baseline', baselinemenuopt)

        procmenuopt = {'Process Data':self.process_data,
                       'Peak Fit': self.peak_fit}
        self.make_menu(menubar, 'Process', procmenuopt)

        revmenuopt = {'Review fitting':self.review}
        self.make_menu(menubar, 'Review', revmenuopt)

        plotmenuopt = {'Plot single spectra':self.file_open,
                       'Plot line data':self.plot_ls,
                       'Plot map data':self.ask_map,
                       'Plot histograms': self.plot_histogram}
        self.make_menu(menubar, 'Plot', plotmenuopt)
        
        manualmenuopt = {'Fit N region manually': self.man_N_fit,
                         'Fit peak manually': self.man_peakfit}
        self.make_menu(menubar, 'Manual fit', manualmenuopt)
        
        manualmenuopt = {'Model N aggregation': self.N_2stage}
        self.make_menu(menubar, '2-stage modelling', manualmenuopt)

        helpmenuopt = {'Read Manual':(lambda s=self: webbrowser.open(s.github_url)),
                       'About':self.about}
        helpmenu = self.make_menu(menubar, 'Help', helpmenuopt)
        helpmenu.insert_separator(1)

        self.master.config(menu=menubar)


##############################################################################

        row = 0

        self.main_fig = Figure(dpi=100)
        self.ax = self.main_fig.add_subplot(111)
        self.main_fig.gca().invert_xaxis()
        self.main_fig.suptitle('QUIDDIT')
        self.IIa_spec = np.loadtxt(settings.IIa_path, delimiter=',')
        self.ax.plot(self.IIa_spec[:, 0], self.IIa_spec[:, 1], 'k-')

        self.main_canvas = self.make_mplcanvas(self.main_fig, erow=row, ecol=0, cspan=4)

        self.message = tk.Text(self, state='disabled', relief=STDRELIEF)
        self.message.grid(row=row, column=4, columnspan=2, sticky=tk.NSEW, padx=5)
        scrl_bar = tk.Scrollbar(self, command=self.message.yview)
        scrl_bar.grid(row=row, column=6, sticky=tk.NSEW)
        self.message['yscrollcommand'] = scrl_bar.set

        row += 1

        toolbar_frame = tk.Frame(self)
        toolbar_frame.grid(row=row, column=1, sticky=tk.W)
        self.toolbar = NavigationToolbar2TkAgg(self.main_canvas, toolbar_frame)
        #self.toolbar = CustomToolbar(self.main_canvas, toolbar_frame)

        self.position = 0

        self.Back_button = self.makebutton(erow=row, ecol=0, caption='Previous',
                                           cmd=lambda s=self: s.plot_seq('prev'),
                                           state='disabled', padx=5, pady=5)
        self.Next_button = self.makebutton(erow=row, ecol=3, caption='Next',
                                           cmd=lambda s=self: s.plot_seq('next'),
                                           state='disabled', padx=5, pady=5, sticky='e')

        self.Histo_button = self.makebutton(erow=row, ecol=2, caption='Histogram',
                                           cmd=self.histo_frame,
                                           state='disabled', padx=5, pady=5, sticky='e')

        for i in range(row):
            self.rowconfigure(i, weight=1, pad=5)
        self.rowconfigure(row, weight=0, pad=5)

        for j in range(6):
            self.columnconfigure(j, weight=1, pad=5)

        self.pack(fill=tk.BOTH, expand=tk.YES)


###############################################################################
    def about(self):
        """Create a toplevel window with version and contact information
        """
        self.toplevel = tk.Toplevel()
        self.toplevel.title('About')

        self.logo_img = tk.PhotoImage(file=r'C:/Users\ls13943\Dropbox\coding\QUIDDIT\QUIDDIT logo.gif')
        self.logo = tk.Label(self.toplevel, image=self.logo_img)
        self.logo.grid(row=0, column=0, padx=5, pady=5)

        about_msg = 'version {}\n(Laura Speich, 10/2017)\n\nlaura.speich@bristol.ac.uk'.format(str(QUIDDITversion))
        msg = tk.Label(self.toplevel, text=about_msg)
        msg.grid(row=1, column=0, padx=5, pady=5)
        Dismiss_button = tk.Button(self.toplevel, text='Dismiss',
                                   command=self.toplevel.destroy,
                                   height=1, width=6, default='active')
        Dismiss_button.grid(row=2, column=0, padx=5, pady=5)
        
        #this doesn't seem to work:
        #for i in range(3):
        #    self.toplevel.rowconfigure(i, weight=1, pad=5)
        #self.toplevel.columnconfigure(1, weight=1, pad=5)
        
        #self.toplevel.grid_columnconfigure(0, weight=1)
        #self.toplevel.grid_rowconfigure(0, weight=1)
        #self.toplevel.resizable(True, True)


    def baseline(self):
        """Baseline correct data
        """
        self.selected_files = fd.askopenfilenames(parent=self, initialdir=self.home,
                                                  title='Select spectra (CSV files)',
                                                  filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        self.file_count.set('{} file(s) selected'.format(len(self.selected_files)))
        if self.selected_files:
            self.save_dir = tk.filedialog.askdirectory(parent=self, initialdir=self.home,
                                                       title='Select directory to save corrected files')
            if self.save_dir:
                loading = LoadingFrame(self.master, len(self.selected_files))
                for item in self.selected_files:
                    self.print_message(self.message,
                                       'Processing file no. {}: {}\n'.format(str(loading.progress['value']+1), str(item)))
                    loading.progress['value'] += 1
                    self.clear_plot(self.main_fig)
                    self.ax = self.main_fig.add_subplot(111)
                    self.plot_spec(item, self.main_fig, self.ax)
                    self.main_canvas.draw()
                    self.update()
                    QUIDDIT_baseline.main(item, self.save_dir)
                loading.destroy()
                self.print_message(self.message, 'Done.')


    def clear_plot(self, fig):
        """Clear all axes in fig
        """
        for ax in fig.get_axes():
            fig.delaxes(ax)
        fig.suptitle('QUIDDIT')
        

    def click_cancel(self):
        """
        """
        self.restore_Ndefault()
        self.N_comp = self.Cvar.get(), self.Avar.get(), self.Xvar.get(), self.Bvar.get(), self.Dvar.get(), self.constvar.get()
        self.toplevel.destroy()


    def click_exit(self):
        """
        """
        #print('The user clicked "Exit"')
        self.master.destroy()


    def click_ok(self):
        """
        """
        print('The user clicked "Ok"')


    def click_save(self):
        """
        """
        self.N_comp = self.Cvar.get(), self.Avar.get(), self.Xvar.get(), self.Bvar.get(), self.Dvar.get(), self.constvar.get()
        self.toplevel.destroy()

    
    def convert_ENVI(self):
        self.selected_ENVI = fd.askopenfilename(parent=self, initialdir=self.home,
                                                title='Select ENVI file (.dat) to open')
        path = os.path.dirname(self.selected_ENVI)
        self.selected_hdr = fd.askopenfilename(parent=self, initialdir=path,
                                               title='Select header (.hdr) file to open')
        #self.data_spacing = sd.askfloat('Spatial resolution', 
        #                                             prompt='What is the spacing between data points (microns)?')
        self.target_directory = fd.askdirectory(parent=self, initialdir=path,
                                                title='Select directory for individual spectra to be stored.')

        #if (self.selected_ENVI and self.selected_hdr and self.data_spacing and self.target_directory):
        if (self.selected_ENVI and self.selected_hdr and self.target_directory):
            #hdr = open(self.selected_hdr, 'r+')
            #hdr.seek(0)
            
            with open(self.selected_hdr, 'r+') as hdr:
                bo_found = any('byte order' in line for line in hdr)
            if not bo_found:
                self.print_message(self.message,
                                       'Byte order not found.')
                bo = self.ask_byteorder()
                self.print_message(self.message,
                                       'Adding byte order = {} to file.'.format(bo))

                hdr = open(self.selected_hdr, 'a+')
                #hdr.seek(0, os.SEEK_END)
                hdr.write('byte order = {}\n'.format(bo))
                hdr.close()
            
            
            
            #for line in hdr.readlines():
            #    if 'byte order' in line:
            #        break
            #    else:
                    #if sys.byteorder == 'little':
                    #    bo = 0
                    #else:
                    #    bo = 1
                    #self.print_message( self.message,
                                       #'Byte order not found.\nAdding byte = {} order to header file'.format(bo))
             #       self.print_message( self.message,
             #                          'Byte order not found.')
              #      bo = self.ask_byteorder()
                    
                    
               #     hdr = open(self.selected_hdr, 'a+')
                #    hdr.write('byte order = {}\n'.format(bo))
                 #   hdr.close()
            
            envi_img = envi.open(self.selected_hdr, self.selected_ENVI).load()
            img_data = envi_img.load()
            xspacing = np.asarray(envi.read_envi_header(self.selected_hdr)['pixel size'], dtype=float)[0]
            yspacing = np.asarray(envi.read_envi_header(self.selected_hdr)['pixel size'], dtype=float)[1]
            wavenum = np.asarray(envi.read_envi_header(self.selected_hdr)['wavelength'], dtype=float)
        
            rows = np.shape(img_data)[0]
            columns = np.shape(img_data)[1]
            print(rows*columns)
        
            loading = LoadingFrame(self.master, (rows*columns))
            loading.progress['value'] = 0
            self.print_message(self.message, 'converting ENVI to CSV. Files will be stored here: {}'.format(self.target_directory))
            
            for i in range(rows):
                for j in range(columns): 
                    spectrum = np.column_stack((wavenum, img_data[i,j,:].flatten()))
                    self.print_message(self.message, 'Saving spectrum {} of {}.'.format(loading.progress['value']+1, rows*columns))
                    #x = i*self.data_spacing
                    #y = j*self.data_spacing
                    x = i * xspacing
                    y = j * yspacing
                    fname = 'X{} Y{}.CSV'.format(str(x), str(y))
                    np.savetxt((self.target_directory + '/' + fname), spectrum, delimiter=',')
                    loading.progress['value'] += 1
                    
                    self.update()

            loading.destroy()
            self.print_message(self.message, 'Done.')


    def file_open(self):
        """Open and plot spectrum
        """
        self.selected_items = fd.askopenfilenames(parent=self, initialdir=self.home,
                                                  title='Select spectra (CSV) to display',
                                                  filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        self.print_message(self.message,
                           'Opening the following files:\n')
        for item in self.selected_items:
            self.print_message(self.message,
                               item)

        self.file_count.set('{} file(s) selected'.format(len(self.selected_items)))

        if self.selected_items:
            self.plotmode.set('single')
            self.plot_seq('zero')


    def goto_github(self):
        """Open default browser and go to github repo
        """
        webbrowser.open(self.github_url)


    def hello(self):
        """for testing
        """
        self.print_message(self.message, "Sorry, this doesn't do anything yet")


    def histo_frame(self):
        row = 0
        self.toplevel = QUIDDITToplevel('Histogram')

        self.histo_fig = Figure(dpi=100)
        self.histo_ax = self.histo_fig.add_subplot(111)
        self.histo_ax.hist(self.plot_item[~np.isnan(self.plot_item)], bins=100)
        self.histo_fig.suptitle(self.map_title)
        self.histo_canvas = self.toplevel.make_mplcanvas(self.histo_fig, erow=row, ecol=0,
                                                   rspan=3)

        tk.Label(self.toplevel, text='Change min/max values for map').grid(row=row, padx=5, pady=5, column=1, columnspan=2)

        row += 1
        self.min = self.toplevel.make_double_entry(lcol=1, lrow=row, ecol=2, erow=row,
                                          caption='min',
                                          textvariable=self.minvar)

        row += 1
        self.max = self.toplevel.make_double_entry(lcol=1, lrow=row, ecol=2, erow=row,
                                          caption='max',
                                          textvariable = self.maxvar)

        row += 1
        self.toplevel.makebutton(erow=row, ecol=1, cspan=2,
                                 padx=5, pady=5, sticky=tk.NSEW,
                                 caption='Redo map',
                                 cmd=lambda s=self: s.redo_map(self.map_grid,
                                                               self.map_title,
                                                               self.histo_canvas,
                                                               self.histo_fig,
                                                               self.histo_ax,
                                                               self.extent))

        toolbar_frame = tk.Frame(self.toplevel)
        toolbar_frame.grid(row=row, column=0)
        self.toolbar = NavigationToolbar2TkAgg(self.histo_canvas, toolbar_frame)
        

    def input_frame(self):
        """Create and input frame to retrieve user input
        """
        self.toplevel = QUIDDITToplevel('Input')
        self.toplevel.bind('<Return>', self.toplevel.destroy)

        row = 0
        tk.Label(self.toplevel,
                 text="Please enter sample name and mantle storage duration").grid(row=row, column=0, columnspan=3, sticky='w')

        row += 1
        self.sample_name = self.toplevel.makeentry(lrow=row, erow=row,
                                                   caption="Sample name",
                                                   width=24,
                                                   textvariable=self.namevar)
        self.sample_name.bind('<Tab>', self.on_input)
        self.sample_name.focus_force()

        row += 1
        self.result_name = self.toplevel.makeentry(erow=row, lrow=row,
                                                   caption='Name for results file: ',
                                                   width=24,
                                                   textvariable=self.resultvar)
        self.set_entry_text(self.result_name, '[sample name] results')
        tk.Label(self.toplevel, text='.csv').grid(row=row, column=2)

        row += 1
        self.review_name = self.toplevel.makeentry(erow=row, lrow=row,
                                                   caption='Name for review file: ',
                                                   width=24,
                                                   textvariable=self.reviewvar)
        self.set_entry_text(self.review_name, '[sample name] review')
        tk.Label(self.toplevel, text='.csv').grid(row=row, column=2)

        row += 1
        self.age = self.toplevel.make_double_entry(erow=row, lrow=row,
                                                   caption='storage duration',
                                                   width=24,
                                                   textvariable=self.agevar)

        tk.Label(self.toplevel, text='(Ma)').grid(row=row, column=2)

        row += 1
        tk.Label(self.toplevel, text='hint: use <tab> to auto-complete file names').grid(row=row, column=0, columnspan=3, padx=5, pady=5, sticky='nesw')

        row += 1
        self.toplevel.makebutton(erow=row, ecol=1, cspan=3,
                                 width=5,
                                 caption='OK',
                                 cmd=self.toplevel.destroy,
                                 sticky=tk.NSEW)

        self.toplevel.deiconify()
        self.toplevel.wait_window()
        
        if self.reviewvar.get() == self.resultvar.get():
            self.reviewvar = self.getvar(self.reviewvar.get()+'2')
            self.print_message(self.message, 'Warning: result and review file have the same name. Saving review as {}.'.format(self.reviewvar.get()))
        
        self.user_inp = (self.namevar.get(), self.resultvar.get(),
                         self.reviewvar.get(), self.agevar.get())

        return self.user_inp
    
    def twostinput_frame(self):
        """Create and input frame to retrieve user input
        """
        self.toplevel = QUIDDITToplevel('2-stage modelling')
        self.toplevel.bind('<Return>', self.toplevel.destroy)

        row = 0
        tk.Label(self.toplevel,
                 text='Please enter data for core and rim').grid(row=row, column=0, columnspan=3, sticky='w')

        row += 1

        self.age = self.toplevel.makeentry(lrow=row, erow=row,
                                                   caption='total duration',
                                                   width=24,
                                                   textvariable=self.agevar)
        tk.Label(self.toplevel, text='(Ma)').grid(row=row, column=2)

        row += 1
        
        tk.Label(self.toplevel,
                 text='core:').grid(row=row, column=0, columnspan=3, sticky='w')
        
        row += 1

        self.c_NT = self.toplevel.makeentry(erow=row, lrow=row,
                                                   caption='[NT]: ',
                                                   width=24,
                                                   textvariable=self.c_NT_var)
        tk.Label(self.toplevel, text='ppm').grid(row=row, column=2)

        row += 1

        self.c_agg = self.toplevel.makeentry(erow=row, lrow=row,
                                                   caption='prop. of B',
                                                   width=24,
                                                   textvariable=self.c_agg_var)
        tk.Label(self.toplevel, text='[-]').grid(row=row, column=2)

        row += 1
        
        tk.Label(self.toplevel,
                 text='rim:').grid(row=row, column=0, columnspan=3, sticky='w')
        
        row += 1

        self.r_NT = self.toplevel.makeentry(erow=row, lrow=row,
                                                   caption='[NT]: ',
                                                   width=24,
                                                   textvariable=self.r_NT_var)
        tk.Label(self.toplevel, text='ppm').grid(row=row, column=2)

        row += 1

        self.r_agg = self.toplevel.makeentry(erow=row, lrow=row,
                                                   caption='prop. of B',
                                                   width=24,
                                                   textvariable=self.r_agg_var)
        tk.Label(self.toplevel, text='[-]').grid(row=row, column=2)

        row += 1

        row += 1
        self.toplevel.makebutton(erow=row, ecol=1, cspan=3,
                                 width=5,
                                 caption='OK',
                                 cmd=self.toplevel.destroy,
                                 sticky=tk.NSEW)

        self.toplevel.deiconify()
        self.toplevel.wait_window()
        user_inp = (self.agevar.get(),
                    self.c_NT_var.get(), self.c_agg_var.get(),
                    self.r_NT_var.get(), self.r_agg_var.get())

        return user_inp

    def peakfit_inp(self):
        self.toplevel = QUIDDITToplevel('Peak Fit')
        self.toplevel.bind('<Return>', self.toplevel.destroy)

        row = 0
        tk.Label(self.toplevel,
                 text='Please enter data for peak fitting').grid(row=row, column=0, columnspan=3, sticky='w')
        
        row +=1
        self.sample_name = self.toplevel.makeentry(lrow=row, erow=row,
                                                   caption="Sample name",
                                                   width=24,
                                                   textvariable=self.namevar)
        self.sample_name.bind('<Tab>', self.on_input2)
        self.sample_name.focus_force()

        row += 1 
        self.selected_peak = self.toplevel.makeentry(lrow=row, erow=row,
                                                   caption='approx. wavenumber',
                                                   width=24,
                                                   textvariable=self.peakvar)
        tk.Label(self.toplevel, text='(cm-1)').grid(row=row, column=2)
        self.selected_peak.bind('<Tab>', self.on_input2)

        row += 1
        self.result_name = self.toplevel.makeentry(erow=row, lrow=row,
                                                   caption='Name for results file: ',
                                                   width=24,
                                                   textvariable=self.resultvar)
        
        
        self.set_entry_text(self.result_name, '[sample name] peak fit')
        tk.Label(self.toplevel, text='.csv').grid(row=row, column=2)
        
        row += 1
        self.toplevel.makebutton(erow=row, ecol=1, cspan=3,
                                 width=5,
                                 caption='OK',
                                 cmd=self.toplevel.destroy,
                                 sticky=tk.NSEW)
        
        self.toplevel.deiconify()
        self.toplevel.wait_window()
        user_inp = (self.namevar.get(),
                    self.peakvar.get(),
                    self.resultvar.get())
        
        return user_inp
        
        

    def loaded(self):
        """Das Window wurde aufgebaut
        """
        self.print_message(self.message, 'Welcome to QUIDDIT ver. 2.0\n')
        self.std = settings.std

    def man_N_fit(self):        
        self.selected_items = fd.askopenfilename(parent=self, initialdir=self.home,
                                                 title='Select spectrum (CSV) files to process',
                                                 filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        
        self.age = tk.simpledialog.askfloat('Manual N aggregation fit', prompt='What is the mantle storage duration (Ma)?', initialvalue=2900.)
        
        if self.selected_items:
            self.clear_plot(self.main_fig)
            spectrum = np.loadtxt(self.selected_items, delimiter=',')
            self.ax = self.main_fig.add_subplot(111)
            self.ax.invert_xaxis()
            self.main_fig.suptitle(self.selected_items.split('/')[-1])

            fit_area = utility.spectrum_slice(spectrum, 1001, 1399)
            self.wav_new = np.arange(fit_area[0][0], fit_area[-1][0], 0.01)
            self.fit_area_inter = utility.inter(fit_area, self.wav_new)
            self.fit_area_inter = self.fit_area_inter.flatten()
            
            C = np.column_stack((settings.std[:,0], settings.std[:,1]))
            A = np.column_stack((settings.std[:,0], settings.std[:,2]))    #generate C, A, X, B and D std
            X = np.column_stack((settings.std[:,0], settings.std[:,3]))    #spectra from CAXBD file
            B = np.column_stack((settings.std[:,0], settings.std[:,4]))   
            D = np.column_stack((settings.std[:,0], settings.std[:,5]))

            C_new = utility.inter(C, self.wav_new)   
            A_new = utility.inter(A, self.wav_new)
            X_new = utility.inter(X, self.wav_new)                     # interpolate C, A, X, B and D spectra
            B_new = utility.inter(B, self.wav_new)
            D_new = utility.inter(D, self.wav_new)
            self.all_comp = np.column_stack((C_new, A_new, X_new, B_new, D_new))
                        
            if self.N_comp[-1] == 1:
                polyx0 = fit_area[-1,1]    
                if polyx0 >0:        
                    polybounds = (0., polyx0)    
                else:
                    polybounds = (polyx0, 0.)
            else:
                polyx0 = 0
                polybounds = (0.,0.)
            
            x0 = [i for i,j in zip((.5, .5, .1, .5, 0., -polyx0), self.N_comp) if j==1]
            bounds =  [i for i,j in zip([(0.,None),(0.,None),(0.,None),(0.,None),(0., None), polybounds], self.N_comp) if j==1]
            
            fit_args = self.all_comp[:,np.where(self.N_comp[:-1]==1)[0]]
            
            fit_res = op.minimize(utility.CAXBD_err, x0=x0, args=(fit_args, self.fit_area_inter), method='SLSQP', bounds=bounds)
            print(fit_res)
            fit = utility.CAXBD(fit_res.x, fit_args)

            self.main_fig.subplots_adjust(left=0.25, bottom=0.4)

            self.ax.plot(fit_area[:,0], fit_area[:,1], 'k.', label='data')
            self.l, = self.ax.plot(self.wav_new, fit, 'g-')
            self.l2, = self.ax.plot(self.wav_new, (fit - self.fit_area_inter),'r-')
            self.ax.axhline(y=0, ls='--', color='k')

            axcolor = 'lightgoldenrodyellow'
            
            ax_C = self.main_fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
            ax_A = self.main_fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
            ax_X = self.main_fig.add_axes([0.25, 0.2, 0.65, 0.03], axisbg=axcolor)
            ax_B = self.main_fig.add_axes([0.25, 0.25, 0.65, 0.03], axisbg=axcolor)
            ax_D = self.main_fig.add_axes([0.25, 0.3, 0.65, 0.03], axisbg=axcolor)
            ax_poly1 = self.main_fig.add_axes([0.25, 0.35, 0.65, 0.03], axisbg=axcolor)
            #axes = [i for i,j in zip((ax_C, ax_A, ax_X, ax_B, ax_D, ax_poly1), self.N_comp) if j==1]
            axes = (ax_C, ax_A, ax_X, ax_B, ax_D, ax_poly1)

            self.sliders = []

            
            #names = [i for i,j in zip(('C', 'A', 'X', 'B', 'D', 'const.'), self.N_comp) if j==1]
            limits = {'C':(0,10), 'A':(0,20), 'X':(0,2), 'B':(0,3), 'D':(0,3), 'const.':(1,1)}           
            
            
            i = 0
                
            for yesno, axis, name in zip(self.N_comp, axes, limits):
                if name == 'const.':
                    slider = Slider(axis, name, fit_res.x[3]-5, fit_res.x[3]+5, valinit = fit_res.x[3], valfmt='%1.1f')
                else:
                    if yesno == 1:
                        slider = Slider(axis, name, limits[name][0]*fit_res.x[i], limits[name][1]*fit_res.x[i], valinit=fit_res.x[i])
                        i += 1
                    else:
                        slider = Slider(axis, name, 0, 5, valinit=0, valfmt='%1.1f')
                slider.on_changed(self.N_widget_update)
                self.sliders.append(slider)
            
            C = self.sliders[0].val
            A = self.sliders[1].val
            X = self.sliders[2].val
            B = self.sliders[3].val
            D = self.sliders[4].val
            
            N_c = np.round(C * 25, 1)
            N_a = np.round(A * 16.5, 1)
            N_b = np.round(B * 79.4, 1)
            N_t = N_a + N_b + N_c
            IaB = np.round(N_b/N_t)

            self.T = np.round(utility.Temp_N(self.age*1e6*365*24*60*60, N_t, N_b/N_t))
            self.fig_text = self.main_fig.text(0.01, 0.75,
                                               '[NC]: {}\n[NA]: {}\n[NB]: {}\n%B.: {}\nT: {}C \nmax D: {}'.format(N_c, N_a, N_b, IaB, self.T, np.round(B*0.365, 2)))             

            self.resetax = self.main_fig.add_axes([0.8, 0.025, 0.1, 0.04])

            self.reset_button = mplButton(self.resetax, 'Reset', color=axcolor, hovercolor='0.975')
            self.reset_button.on_clicked( lambda event, arg=self.sliders: self.widget_reset(event, arg))
          
            
        

    def man_peakfit(self):
        
        self.selected_items = fd.askopenfilename(parent=self, initialdir=self.home,
                                                  title='Select spectrum (CSV) to process',
                                                  filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        self.peak = tk.simpledialog.askfloat('Manual peak fit', prompt='Where is your peak located (approx. wavenumber in cm-1)?')
        
        if (self.selected_items and self.peak !=0):
            self.clear_plot(self.main_fig)
            lower_bound = self.peak - 30
            upper_bound = self.peak + 30
            
            spectrum = np.loadtxt(self.selected_items, delimiter=',')
            self.ax = self.main_fig.add_subplot(111)
            self.main_fig.suptitle(self.selected_items.split('/')[-1])

            fit_area = utility.spectrum_slice(spectrum, lower_bound, upper_bound)
            self.wav_new = np.arange(fit_area[0][0], fit_area[-1][0], 0.01)
            self.fit_area_inter = utility.inter(fit_area, self.wav_new)

            pos_guess = fit_area[:,0][np.argmax(fit_area[:,1])]
            height_guess = fit_area[:,1][np.argmax(fit_area[:,1])]
            width_guess_l = 2
            width_guess_r = 2
            sigma_guess = 1
            const_guess = fit_area[-1,1]

            #x0=[(pos_guess, height_guess*1.2, width_guess_l, width_guess_r, sigma_guess)]
            x0=[(pos_guess, height_guess*1.2, width_guess_l, width_guess_r, sigma_guess, const_guess)]
            bounds = [(pos_guess-3,pos_guess+3),(0.0,None),(0.0, None),(0.0,None), (0,1), (None,None)]
            fit_args = (self.wav_new, self.fit_area_inter)

            
            # optimization: 
            #fit_res = op.minimize(utility.pseudovoigt, x0=x0, args=fit_args, method='SLSQP', bounds=bounds)
            fit_res = op.minimize(utility.pseudovoigt_const, x0=x0, args=fit_args, method='SLSQP', bounds=bounds)

            
            #fit = utility.pseudovoigt_fit(self.wav_new, *fit_res.x)
            fit = utility.pseudovoigt_fit(self.wav_new, *fit_res.x[:-1]) + fit_res.x[-1]


            self.main_fig.subplots_adjust(left=0.25, bottom=0.4)


            self.ax.plot(fit_area[:,0], fit_area[:,1], 'k.', label='data')
            self.ax.invert_xaxis()
            self.l, = self.ax.plot(self.wav_new, fit, 'g-')
            self.l2, = self.ax.plot(self.wav_new, (fit - self.fit_area_inter),'r-')
            self.ax.axhline(y=0, ls='--', color='k')

            axcolor = 'lightgoldenrodyellow'
            
            ax_const = self.main_fig.add_axes([0.25, 0.35, 0.65, 0.03], facecolor=axcolor)
            ax_x0 = self.main_fig.add_axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)
            ax_I = self.main_fig.add_axes([0.25, 0.25, 0.65, 0.03], facecolor=axcolor)
            ax_HWHM_l = self.main_fig.add_axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
            ax_HWHM_r = self.main_fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
            ax_sigma = self.main_fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
            

            self.s_x0 = Slider(ax_x0, 'peak pos.', pos_guess-3, pos_guess+3, valinit=fit_res.x[0], valfmt='%1.1f')
            self.s_I = Slider(ax_I, 'peak height', 0, fit_res.x[1]+fit_res.x[1]*0.25, valinit=fit_res.x[1], valfmt='%1.1f')
            self.s_HWHM_l = Slider(ax_HWHM_l, 'l. half width', 0, fit_res.x[2]*3, valinit=fit_res.x[2], valfmt='%1.1f')
            self.s_HWHM_r = Slider(ax_HWHM_r, 'r. half width', 0, fit_res.x[3]*3, valinit=fit_res.x[3], valfmt='%1.1f')
            self.s_sigma = Slider(ax_sigma, 'Lorentz. contr.', 0, 1, valinit = fit_res.x[4], valfmt='%1.1f')
            self.s_const = Slider(ax_const, 'const.', fit_res.x[-1]*0.7, fit_res.x[-1]*1.3, valinit = fit_res.x[-1], valfmt='%1.1f')


            self.fig_text = self.main_fig.text(0.29, 0.8,
                                  'Peak area:\n{} cm-2'.format(np.round(utility.peak_area(self.s_I.val,
                                               self.s_HWHM_l.val, self.s_HWHM_r.val,
                                               self.s_sigma.val)), 2))

            self.s_x0.on_changed(self.p_widget_update)
            self.s_I.on_changed(self.p_widget_update)
            self.s_HWHM_l.on_changed(self.p_widget_update)
            self.s_HWHM_r.on_changed(self.p_widget_update)
            self.s_sigma.on_changed(self.p_widget_update)
            self.s_const.on_changed(self.p_widget_update)
            self.resetax = self.main_fig.add_axes([0.8, 0.025, 0.1, 0.04])

            self.reset_button = mplButton(self.resetax, 'Reset', color=axcolor, hovercolor='0.975')

            #sliders = (self.s_x0, self.s_I, self.s_HWHM_l, self.s_HWHM_r, self.s_sigma)
            sliders = (self.s_x0, self.s_I, self.s_HWHM_l, self.s_HWHM_r, self.s_sigma, self.s_const)
            self.reset_button.on_clicked( lambda event, arg=sliders: self.widget_reset(event, arg))

    def N_2stage(self):
        self.user_inp = self.twostinput_frame()
        QUIDDIT_two_stage_model.main(self.user_inp[0], self.user_inp[1],
                                     self.user_inp[2], self.user_inp[3],
                                     self.user_inp[4])
        
        self.clear_plot(self.main_fig)
        self.plot_2st('output.csv', self.main_fig)
        self.print_message(self.message, '2-stage modelling:')
        self.print_message(self.message, 'total duration: %s Ma\ncore [N]t: %s ppm\ncore agg.: %s\nrim [N]t: %s ppm\nrim agg.: %s' %self.user_inp)

    def on_input(self, event):
        """Get sample name and use it to suggest file names
        """
        inp = self.namevar.get()
        self.result_name.delete('0', 'end')
        self.result_name.insert('0', inp+' results')
        self.review_name.delete('0', 'end')
        self.review_name.insert('0', inp+' review')
        
    def on_input2(self, event):
        """Get sample name and use it to suggest file names
        """
        name = self.namevar.get()
        pk = str(self.peakvar.get()).replace('.','pt')
        self.result_name.delete('0', 'end')
        self.result_name.insert('0', '{} {} peak fit'. format(name, pk))


    def plot_ls(self):
        """Method to plot results of linescans.
        """
        self.res_file = fd.askopenfilename(parent=self, initialdir=self.home,
                                           title='Select results file',
                                           filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        self.clear_plot(self.main_fig)
        self.main_fig.suptitle('')

        res = np.loadtxt(self.res_file, delimiter=',', dtype=utility.results_dtype, skiprows=2)
        points = []
        for item in res['name']:
            points.append(float(item.split()[-1][-8:-4]))

        self.ax0 = self.main_fig.add_subplot(3, 3, 1)
        self.ax0.text(0, 0.5, self.res_file.split('/')[-1], bbox={'facecolor':'white', 'pad':10})
        self.ax0.axis('off')

        self.ax1 = self.main_fig.add_subplot(3, 3, 2)
        self.ax1.plot(res['[NB]']/79.4, res['p_area_ana'], 'k.')
        self.ax1.plot((0, 100), (0, 6400), 'k-')
        self.ax1.set(xlim=(0, 10), ylim=(0, 640),
                     xlabel='$\mathregular{\mu_B}$', ylabel="I(B') $\mathregular{(cm^{-2})}$")

        self.ax2 = self.main_fig.add_subplot(3, 3, 3)
        self.ax2.plot(points, res['H_area_ana'], 'k.')
        self.ax2.set(ylim=(0, None),
                     ylabel='I(3107) $\mathregular{(cm^{-2})}$',
                     xticklabels=[])

        self.ax3 = self.main_fig.add_subplot(3, 3, 4)
        self.ax3.plot(points, res['[NA]'], '.', label='$\mathregular{[N_A]}$')
        self.ax3.plot(points, res['[NB]'], '.', label='$\mathregular{[N_B]}$')
        self.ax3.plot(points, res['[NT]'], '.', label='$\mathregular{[N_t]}$')
        self.ax3.set(ylim=(0, None),
                     ylabel='concentration (ppm)',
                     xticklabels=[])
        self.ax3.legend(loc='best')

        self.ax4 = self.main_fig.add_subplot(3, 3, 5)
        self.ax4.plot(points, (res['[NB]']/(res['[NA]']+res['[NB]'])), 'k.')
        self.ax4.set(ylim=(0, 1),
                     ylabel='$\mathregular{[N_B]/[N_T]}$',
                     xticklabels=[])

        self.ax5 = self.main_fig.add_subplot(3, 3, 6)
        self.ax5.plot(points, res['T'], 'k.')
        self.ax5.set(ylim=(1050, 1350),
                     ylabel='$\mathregular{T_N (^{\circ}C)}$',
                     xticklabels=[])

        self.ax6 = self.main_fig.add_subplot(3, 3, 7)
        self.ax6.plot(points, res['p_area_ana'], 'k.')
        self.ax6.set(ylim=(0, None),
                     xlabel='#', ylabel="I(B') $\mathregular{(cm^{-2})}$")

        self.ax7 = self.main_fig.add_subplot(3, 3, 8)
        P0 = (res['[NB]']/79.4)*64
        deg = (1-(res['p_area_ana']/P0))
        self.ax7.plot(points, deg, 'k.')
        self.ax7.set(ylim=(0, 1),
                     xlabel='#', ylabel="$\mathregular{1-(I(B')/I(B')_0)}$")

        self.ax8 = self.main_fig.add_subplot(3, 3, 9)
        lnstuff = np.log(np.log(P0/res['p_area_ana'])/(1000*1e6*365*24*60*60))
        TP = ((88446)/(19.687-lnstuff)) -273
        self.ax8.plot(points, TP, 'k.')
        self.ax8.set(ylim=(1050, 1350),
                     xlabel='#', ylabel='$\mathregular{T_P (^{\circ}C)}$')

        self.main_fig.tight_layout(pad=0.7, w_pad=0, h_pad=0)
        self.main_canvas.draw()

    def ask_byteorder(self):
        self.toplevel = tk.Toplevel()
        self.toplevel.title('Error. Byte order not found.')
        tk.Label(self.toplevel, text='Was the ENVI file created on Windows or MacOS/Linux?').grid(row=0, column=0, columnspan=2, padx=5, pady=5)
        byteorder = self.getvar(0)
        tk.Radiobutton(self.toplevel, text='Windows', variable=byteorder, value=0).grid(row=1, column=0, padx=5, pady=5)
        tk.Radiobutton(self.toplevel, text='MacOS/Linux', variable=byteorder, value=1).grid(row=1, column=1, padx=5, pady=5)
        btn = tk.Button(self.toplevel, text='Ok', command=self.toplevel.destroy)
        btn.grid(row=2, column=0, padx=5, pady=5)
        self.toplevel.wait_window()
        return byteorder.get()


    def peak_fit(self):
        self.input = self.peakfit_inp()
        
        if self.input:
            (name, peak, filename) = self.input
            resultfile = filename+'.csv'
            self.selected_items = fd.askopenfilenames(parent=self, initialdir=self.home,
                                                      title='Select baseline corrected spectra for peak fitting (CSV)',
                                                      filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
            self.print_message(self.message, 'Fitting peaks at {} for sample {}\nwriting results to {}.\n'.format(peak, name, resultfile))
            
            loading = LoadingFrame(self.master, len(self.selected_items))
            for file in self.selected_items:
                if loading.progress['value'] == 0:
                    with open(resultfile, 'w') as res_fob:
                        res_fob.write('Peak fitting esults for sample {}'.format(name) + ':\n')
                        res_fob.write(utility.peakfit_header + '\n')

                self.clear_plot(self.main_fig)
                self.ax = self.main_fig.add_subplot(111)
                self.plot_spec(file, self.main_fig, self.ax)
                self.main_canvas.draw()
                self.print_message(self.message,
                                   '\nProcessing file no. {} of {}'.format(str(loading.progress['value']+1),
                                                        loading.progress['maximum']))

                peakfit_res = QUIDDIT_peakfit.main(file, peak)
                self.print_message(self.message,
                                   'Results for this spectrum:')
                for index, item in enumerate(zip(peakfit_res, utility.peakfit_header.split(','))):
                    if index == 0:
                        self.print_message(self.message,
                                           '{}: {}'.format(item[1], item[0]))
                    else:
                        self.print_message(self.message,
                                           #'{}: {}'.format(item[1], item[0]))
                                           '{}: {}'.format(item[1], np.round(item[0], 2)))

                with open(resultfile, 'a') as res_fob:
                    for res in peakfit_res:
                        res_fob.write(str(res)+',')
                    res_fob.write('\n')

                loading.progress['value'] += 1
                self.update()

            loading.destroy()
            self.print_message(self.message, 'Done.')
    
    

    def ask_map(self):
        self.res_file = fd.askopenfilename(parent=self, initialdir=self.home,
                                           title='Select results file',
                                           filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        self.selected_items = ('$[N_T]$ (ppm)',
                               '$[N_A]$ (ppm)',
                               '$[N_B]$ (ppm)',
                               '$[N_B]/[N_T]$',
                               '$T (^{\circ}C)$',
                               'platelet peak position $(cm^{-1})$',
                               'platelet peak area $(cm^{-2})$',
                               'platelet peak width $(cm^{-1})$',
                               'platelet peak symmetry $(cm^{-1})$',
                               'I(3107) $(cm^{-2})$')
        self.plotmode.set('map')
        self.plot_seq('zero')

    def plot_histogram(self):
        self.res_file = fd.askopenfilename(parent=self, initialdir=self.home,
                                           title='Select results file',
                                           filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        self.print_message(self.message, 'Opening file: {}'.format(self.res_file))
        self.results = np.loadtxt(self.res_file, delimiter=',', dtype=utility.results_dtype, skiprows=2)

        self.bins = sd.askinteger(
                  "No. of bins",
                  "Enter ingeter value between 3 and 1000",
                  initialvalue=200,
                  minvalue=3,
                  maxvalue=1000)
        self.print_message(self.message, 'Creating histograms with {} bins'.format(self.bins))
        
        self.selected_items = ('$[N_T]$ (ppm)',
                               '$[N_A]$ (ppm)',
                               '$[N_B]$ (ppm)',
                               '$[N_B]/[N_T]$',
                               '$T (^{\circ}C)$',
                               'platelet peak position $(cm^{-1})$',
                               'platelet peak area $(cm^{-2})$',
                               'platelet peak width $(cm^{-1})$',
                               'platelet peak symmetry $(cm^{-1})$',
                               'I(3107) $(cm^{-2})$')
        
        self.plotmode.set('histogram')
        self.plot_seq('zero')


    def plot_map(self, data, title, fig, ax, extent, clim):
        img = ax.imshow(data, origin='lower', extent=extent,
                             cmap=STDCOLS, clim=clim)
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.3)
        self.cbar = fig.colorbar(img, cax=cax)
        fig.suptitle(title)
        return img


    def plot_seq(self, seq):
        """plots spectra within a sequence.
        plotmode=='single': plot single spectrum
        plotmode=='review': plot spectrum and review elements
        plotmode=='map': plot maps
        plotmode=='histogram': plot histograms from results file
        """
        self.clear_plot(self.main_fig)

        self.Next_button['state'] = 'normal'
        self.Back_button['state'] = 'normal'

        if seq == 'zero':
            self.index = 0
        elif seq == 'prev':
            if self.index == 0:
                self.index = len(self.selected_items)-1
            else:
                self.index -= 1
        elif seq == 'next':
            if self.index == len(self.selected_items)-1:
                self.index = 0
            else:
                self.index += 1

        if self.plotmode.get() == 'single':
            self.ax = self.main_fig.add_subplot(1, 1, 1)
            self.plot_spec((self.selected_items[self.index]),
                           self.main_fig, self.ax)

        elif self.plotmode.get() == 'review':
            self.main_fig.suptitle('Review')
            self.print_message(self.message, 'Starting review.')
            for i in range(3):
                self.ax = self.main_fig.add_subplot(3, 1, i+1)
                self.plot_spec(self.selected_items[self.index],
                               self.main_fig, self.ax, xlabel='')
                rev = self.read_review(self.review_path, self.index)
                spec = np.loadtxt(self.selected_items[self.index], delimiter=',')

                N_abs_new = utility.inter(spec, np.array(self.std[:, 0]))
                c = rev['c']
                a = rev['a']
                x = rev['x']
                b = rev['b']
                d = rev['d']
                const = rev['N_poly']
                N_fit = utility.CAXBD(np.array((c, a, x, b, d, const)), self.std[:, 1:])

                NC = np.round(c * 25, 1)
                NA = np.round(a * 16.5, 1)
                NB = np.round(b * 79.4, 1)
                NT = np.sum(np.nan_to_num((NC, NA, NB)))

                p_spec = utility.spectrum_slice(spec, 1327, 1420)
                p_wav = p_spec[:, 0]
                p_params = (rev['p_x0'], rev['p_I'], rev['p_HWHM_l'], rev['p_HWHM_r'], rev['p_sigma'],
                            rev['H1405_x0'], rev['H1405_I'], rev['H1405_HWHM_l'], rev['H1405_HWHM_r'], rev['H1405_sigma'],
                            rev['B_x0'], rev['B_I'], rev['B_HWHM_l'], rev['B_HWHM_r'], rev['B_sigma'],
                            rev['psv_c'])

                p_fit = utility.ultimatepsv_fit(p_wav, *p_params)
                p_peak_area = np.round(utility.peak_area(*p_params[1:5]), 1)

                H_spec = utility.spectrum_slice(spec, 3000, 3200)
                H_wav = H_spec[:, 0]
                H_params = (rev['H_pos'], rev['H_I'], rev['H_HWHM_l'], rev['H_HWHM_r'], rev['H_sigma'])
                H_bg_params = (rev['H_bg_a'], rev['H_bg_b'], rev['H_bg_c'], rev['H_bg_d'])
                H_fit = utility.pseudovoigt_fit(H_wav, *H_params) + np.polyval(H_bg_params, H_wav)
                H_peak_area = np.round(utility.peak_area(*H_params[1:]), 1)

                if i == 0:
                    self.ax.axhline(y=0, color='0.7', linestyle='--')
                    self.ax.plot(np.array(self.std[:, 0]),
                                 N_fit,
                                 'g-', label='fit')
                    self.ax.plot(np.array(self.std[:, 0]),
                                 N_abs_new - N_fit,
                                 'r-', label='misfit')
                    self.ax.set(xlim=(1400, 1000))
                    self.print_message(self.message,
                                       self.selected_items[self.index])
                    self.print_message(self.message,
                                       '\nNC: {} ppm\nNA: {} ppm\nNB: {} ppm\ntotal: {} ppm'.format(NC, NA, NB, NT))
                    self.print_message(self.message,
                                       'platelet peak area: {0:.2f} cm-2'.format(p_peak_area))
                    self.print_message(self.message,
                                       'platelet peak position: {0:.2f} cm-1'.format(rev['p_x0']))
                    self.print_message(self.message,
                                       '3107 peak area: {0:.2f} cm-2\n'.format(H_peak_area))

                elif i == 1:
                    self.ax.axhline(y=0, color='0.7', linestyle='--')
                    self.ax.plot(p_wav, p_fit,
                                 'g-', label='fit')
                    self.ax.plot(p_wav, p_fit - p_spec[:, 1],
                                 'r-', label='misfit')
                    self.ax.set(xlim=(1420, 1327))

                elif i == 2:
                    self.ax.axhline(y=0, color='0.7', linestyle='--')
                    self.ax.plot(H_wav, H_fit,
                                 'g-', label='fit')
                    self.ax.plot(H_wav, H_fit - H_spec[:, 1],
                                 'r-', label='misfit')
                    self.ax.set(xlim=(3200, 3000),
                                xlabel='wavenumber ($\mathregular{cm^{-1}}$)')
                    self.ax.legend(loc='best')

                else:
                    self.print_message(self.message, 'Something went wrong')


        elif self.plotmode.get() == 'map':
            self.print_message(self.message, 'Plotting map. This may take a few seconds...')
            self.Histo_button['state'] = 'normal'
            x = []
            y = []
            res = self.read_results(self.res_file)
            for item in res['name']:
                #xy = item.decode().split('/')[-1]
                xy = item.split('/')[-1]
                x.append(float(xy.split(' ')[0][2:]))
                y.append(float(xy.split(' ')[1].strip(".csv'")[1:]))

            self.extent = (min(x), max(x), min(y), max(y))
            resolution = 2000j


            grid_x, grid_y = np.mgrid[self.extent[0]:self.extent[1]:resolution,
                                      self.extent[2]:self.extent[3]:resolution]

            maps = {'$[N_T]$ (ppm)': res['[NT]'],
                    '$[N_A]$ (ppm)': res['[NA]'],
                    '$[N_B]$ (ppm)': res['[NB]'],
                    '$[N_B]/[N_T]$': (res['[NB]']/res['[NT]']),
                    '$T (^{\circ}C)$': res['T'],
                    'platelet peak position $(cm^{-1})$': res['p_x0'],
                    'platelet peak area $(cm^{-2})$': res['p_area_ana'],
                    'platelet peak width $(cm^{-1})$': (res['p_HWHM_l'] + res['p_HWHM_r']),
                    'platelet peak symmetry $(cm^{-1})$': (res['p_x0'] - res['avg']),
                    'I(3107) $(cm^{-2})$': res['H_area_ana']}


            clims = {'$[N_T]$ (ppm)': (None, None),
                     '$[N_A]$ (ppm)': (None, None),
                     '$[N_B]$ (ppm)': (None, None),
                     '$[N_B]/[N_T]$': (0, 1),
                     '$T (^{\circ}C)$': (1000, 1400),
                     'platelet peak position $(cm^{-1})$':(1358, 1378),
                     'platelet peak area $(cm^{-2})$': (None, None),
                     'platelet peak width $(cm^{-1})$': (None, 25),
                     'platelet peak symmetry $(cm^{-1})$': (-15, 5),
                     'I(3107) $(cm^{-2})$': (None, None)}

            plots = self.selected_items
            clim = clims[plots[self.index]]

            self.plot_item = maps[plots[self.index]]
            self.map_title = plots[self.index]

            if plots[self.index] == '$[N_B]/[N_T]$':
                self.ax.set()

            self.map_grid = utility.make_2dgrid(x, y, grid_x, grid_y, self.plot_item)
            self.ax = self.main_fig.add_subplot(1, 1, 1)
            self.plot_map(self.map_grid, self.map_title, self.main_fig, self.ax, self.extent, clim)
            

        elif self.plotmode.get() == 'histogram':
            
            res = self.read_results(self.res_file)
            histograms = {'$[N_T]$ (ppm)': res['[NT]'],
                    '$[N_A]$ (ppm)': res['[NA]'],
                    '$[N_B]$ (ppm)': res['[NB]'],
                    '$[N_B]/[N_T]$': (res['[NB]']/res['[NT]']),
                    '$T (^{\circ}C)$': res['T'],
                    'platelet peak position $(cm^{-1})$': res['p_x0'],
                    'platelet peak area $(cm^{-2})$': res['p_area_ana'],
                    'platelet peak width $(cm^{-1})$': (res['p_HWHM_l'] + res['p_HWHM_r']),
                    'platelet peak symmetry $(cm^{-1})$': (res['p_x0'] - res['avg']),
                    'I(3107) $(cm^{-2})$': res['H_area_ana']}
            
            plots = self.selected_items
            data = histograms[plots[self.index]]      
            histogram = data[~np.isnan(data)]
            
            self.ax = self.main_fig.add_subplot(1, 1, 1)
            self.ax.hist(histogram, bins=self.bins)
            self.ax.set_title(plots[self.index])

        else:
            self.print_message(self.message, 'Error: Unknown plotmode.')

        self.main_canvas.draw()


    def plot_spec(self, file, fig, ax,
                  xlabel='wavenumber ($\mathregular{cm^{-1}}$)',
                  ylabel='absorption', **options):
        """plot a single spectrum in fig on axis ax by reading a file
        """
        spectrum = np.loadtxt(file, delimiter=',')
        title = file.split('/')[-1]
        ax.plot(spectrum[:, 0], spectrum[:, 1], 'k-', label='data')
        ax.set(xlim=(675, 4000), ylim=(None, None),
               xlabel=xlabel, ylabel=ylabel, **options)
        fig.gca().invert_xaxis()
        fig.suptitle(title)
        

    def plot_2st(self, file, fig):
        data = np.loadtxt(file, delimiter=',', skiprows=1)
        ax = fig.add_subplot(111)
        ax.plot(data[:,0], data[:,1], 'ro', label='core')
        ax.plot(data[:,0], data[:,2], 'bo', label='rim')
        ax.set(ylim=(1100, 1500), xlim=(0, float(self.agevar.get())),
               xlabel='Duration of first anneal (Ma)',
               ylabel ='Temperature ($\mathregular{^{\circ}}$C)')
        ax.legend(loc='best')
        fig.suptitle('2-stage model')
        self.main_canvas.draw()
        

    def print_message(self, textwidget, text):
        """print text to textwidget
        """
        textwidget['state'] = 'normal'
        textwidget.insert('end', text+'\n')
        textwidget['state'] = 'disabled'
        textwidget.see('end')


    def process_data(self):
        """Process data
        """
        self.selected_items = fd.askopenfilenames(parent=self, initialdir=self.home,
                                                  title='Select baseline corrected spectra (CSV)',
                                                  filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        self.input = self.input_frame()

        if self.input:
            self.results = np.zeros((len(self.selected_items),),
                                    dtype=utility.results_dtype)
            self.review = np.zeros((len(self.selected_items),),
                                   dtype=utility.review_dtype)

            self.print_message(self.message,
                               'Starting to process the data.\nsample name: %s\nresults file: %s.csv\nreview file: %s.csv\nsample age: %i Ma' %self.input)
            self.print_message(self.message,
                               'Files with results will be stored here: {}'.format(self.home))

            self.sample, self.res_name, self.rev_name, self.age = self.input

            self.resultfile = self.home + '/' + self.res_name + '.csv'
            self.reviewfile = self.home + '/' + self.rev_name + '.csv'

            i = 0

            loading = LoadingFrame(self.master, len(self.selected_items))
            for item in self.selected_items:
                if loading.progress['value'] == 0:
                    with open(self.resultfile, 'w') as res_fob:
                        res_fob.write('Results for sample %s - age: %.0f Ma' %(str(self.sample), round(self.age, 3)) + ':\n')
                        res_fob.write(utility.res_header+'\n')

                    with open(self.reviewfile, 'w') as rev_fob:
                        rev_fob.write('Review for sample %s' %(str(self.sample)) + ':\n')
                        rev_fob.write(utility.rev_header+'\n')

                self.clear_plot(self.main_fig)
                self.ax = self.main_fig.add_subplot(111)
                self.plot_spec(item, self.main_fig, self.ax)
                self.main_canvas.draw()
                self.print_message(self.message,
                                   '\nProcessing file no. {} of {}'.format(str(loading.progress['value']+1),
                                                        loading.progress['maximum']))

                curr_res, curr_rev = QUIDDIT_main.main(item, self.age, self.N_comp)
                self.print_message(self.message,
                                   'Results for this spectrum:')

                for index, item in enumerate(zip(curr_res[0], utility.res_header.split(','))):
                    if index == 0:
                        self.print_message(self.message,
                                           '{}: {}'.format(item[1], item[0]))
                    else:
                        self.print_message(self.message,
                                           '{}: {}'.format(item[1], np.round(item[0], 2)))

                with open(self.resultfile, 'a') as res_fob:
                    for item in curr_res[0]:
                        res_fob.write(str(item)+',')
                    res_fob.write('\n')

                with open(self.reviewfile, 'a') as rev_fob:
                    for item in curr_rev[0]:
                        rev_fob.write(str(item)+',')
                    rev_fob.write('\n')

                #self.results[i] = curr_res
                #self.review[i] = curr_rev

                loading.progress['value'] += 1
                #loading.update()
                self.update()
                i += 1

            loading.destroy()
            self.print_message(self.message, 'Done.')

        else:
            self.print_message(self.message, 'What do you want me to do?')

        self.set_defaults()
        curr_res, curr_rev = [], []


    def read_results(self, resultfile):
        results = np.loadtxt(resultfile, delimiter=',',
                             dtype=utility.results_dtype, skiprows=2)
        return results


    def read_review(self, reviewfile, seq):
        review = np.loadtxt(reviewfile, delimiter=',',
                            dtype=utility.review_dtype, skiprows=2)
        row = review[seq]
        return row


    def redo_map(self, data, title, canvas, fig, ax, extent):
        fig.delaxes(ax)
        ax = fig.add_subplot(111)
        clim=(self.minvar.get(), self.maxvar.get())
        img = self.plot_map(data, title, fig, ax, extent, clim)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.3)
        self.cbar = fig.colorbar(img, cax=cax)
        canvas.draw()


    def restore_Ndefault(self):
        """Restore default selection of end-members for N fitting
        """
        self.A.select()
        self.B.select()
        self.D.select()
        self.const.select()


    def review(self):
        """
        """
        self.review_path = fd.askopenfilename(parent=self, initialdir=self.home,
                                              title='Select review file',
                                              filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        if self.review_path:
            self.selected_items = fd.askopenfilenames(parent=self, initialdir=self.home,
                                                      title='Select corresponding corrected spectra (CSV)',
                                                      filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
            if self.selected_items:
                self.plotmode.set('review')
                self.plot_seq('zero')


    def set_defaults(self):
        """Set all the defaults
        """
        self.home = settings.home
        self.N_comp = settings.N_comp
        self.file_count = self.getvar('')
        self.github_url = 'https://github.com/LauraSp/QUIDDIT'
        self.namevar = self.getvar('')
        self.resultvar = self.getvar('')
        self.reviewvar = self.getvar('')
        self.agevar = self.getvar(2900)
        self.peakvar = self.getvar(3107.0)
        self.c_NT_var = self.getvar(0.)
        self.r_NT_var = self.getvar(0.)
        self.c_agg_var = self.getvar(0.)
        self.r_agg_var = self.getvar(0.)
        self.plotmode = self.getvar('')
        self.minvar = self.getvar(0.)
        self.maxvar = self.getvar(1.)
        self.peak = self.getvar(0)


    def settings(self):
        """Create a toplevel window to change settings
        """
        self.toplevel = tk.Toplevel()
        self.toplevel.title('Settings')

        blset_frame = tk.LabelFrame(self.toplevel, relief='groove', text='Baseline settings')
        blset_frame.grid(row=0, pady=5, columnspan=2, sticky='w')
        bllabel = tk.Label(blset_frame, text='some settings...')
        bllabel.grid()

        fitset_frame = tk.LabelFrame(self.toplevel, relief='ridge', text='Fit settings')
        fitset_frame.grid(row=1, pady=5, columnspan=2, sticky='w')

        tk.Label(fitset_frame, text='included in N fit:').grid(row=2, column=0, sticky='w')

        self.Cvar = tk.IntVar()

        self.C = tk.Checkbutton(fitset_frame, text='C centre', variable=self.Cvar)
        self.C.grid(sticky='w', row=3, column=0)

        self.Avar = tk.IntVar()
        self.A = tk.Checkbutton(fitset_frame, text='A centre', variable=self.Avar)
        self.A.grid(sticky='w', row=4, column=0)

        self.Xvar = tk.IntVar()
        self.X = tk.Checkbutton(fitset_frame, text='X centre', variable=self.Xvar)
        self.X.grid(sticky='w', row=5, column=0)

        self.Bvar = tk.IntVar()
        self.B = tk.Checkbutton(fitset_frame, text='B centre', variable=self.Bvar)
        self.B.grid(sticky='w', row=3, column=1)

        self.Dvar = tk.IntVar()
        self.D = tk.Checkbutton(fitset_frame, text='D centre', variable=self.Dvar)
        self.D.grid(sticky='w', row=4, column=1)

        self.constvar = tk.IntVar()
        self.const = tk.Checkbutton(fitset_frame, text='add constant', variable=self.constvar)
        self.const.grid(sticky='w', row=5, column=1, pady=3)

        plateletlabel = tk.Label(fitset_frame, text="B' fit")
        plateletlabel.grid(sticky='w', row=6)
        plateletlabel2 = tk.Label(fitset_frame, text='some settings...')
        plateletlabel2.grid(sticky='w', row=7)

        default_button = tk.Button(self.toplevel, text='Restore defaults',
                                   command=self.restore_Ndefault)
        default_button.grid(sticky='nesw', row=8, columnspan=2)
        save_button = tk.Button(self.toplevel, text='Save',
                                command=self.click_save, height=1, width=6, default='active')
        save_button.grid(sticky='w', row=9, column=0, pady=5)
        cancel_button = tk.Button(self.toplevel, text='Cancel',
                                  command=self.click_cancel, height=1, width=6)
        cancel_button.grid(sticky='e', row=9, column=1, pady=5)

        self.restore_Ndefault()

        self.main_fig.canvas.draw_idle()


    def p_widget_update(self, val):
        pos = self.s_x0.val
        I = self.s_I.val
        HWHM_l = self.s_HWHM_l.val
        HWHM_r = self.s_HWHM_r.val
        sigma = self.s_sigma.val
        const = self.s_const.val
        self.l.set_ydata(utility.pseudovoigt_fit(self.wav_new, pos, I, HWHM_l, HWHM_r, sigma )+const)
        self.l2.set_ydata(utility.pseudovoigt_fit(self.wav_new, pos, I, HWHM_l, HWHM_r, sigma )+const-self.fit_area_inter)
        self.fig_text.set(text='Peak area:\n{} cm-2'.format(np.round(utility.peak_area(I,
                                           HWHM_l, HWHM_r, sigma)), 2))
        
    def N_widget_update(self, val):       
        C = self.sliders[0].val
        A = self.sliders[1].val
        X = self.sliders[2].val
        B = self.sliders[3].val
        D = self.sliders[4].val
        poly1 = self.sliders[5].val

        factors = np.array([C, A, X, B, D, poly1])
        
        N_c = np.round(C * 25, 1)
        N_a = np.round(A * 16.5, 1)
        N_b = np.round(B * 79.4, 1)
        N_t = N_a + N_b + N_c
        IaB = np.round(N_b/N_t)
        T = np.round(utility.Temp_N(self.age*1e6*365*24*60*60, N_t, N_b/N_t))
        
        self.l.set_ydata(utility.CAXBD(factors, self.all_comp))
        self.l2.set_ydata(utility.CAXBD(factors, self.all_comp)-self.fit_area_inter)
        self.fig_text.set(text='[NC]: {}\n[NA]: {}\n[NB]: {}\n%B.: {}\nT: {}C \nmax D: {}'.format(N_c, N_a, N_b, IaB, T, np.round(B*0.365, 2)))

        
    def widget_reset(self, event, sliders):
        for slider in sliders:
            slider.reset()


class QUIDDITToplevel(tk.Toplevel, TclWinBase):
    """Making a toplevel window that inherits from TclWinBase
    """
    def __init__(self, title):
        super().__init__()
        self.toplevel = tk.Toplevel()
        self.toplevel.title(title)
        self.toplevel.protocol("WM_DELETE_WINDOW", self.toplevel.destroy)


class LoadingFrame(tk.Frame):
    """Creating a frame for a progress bar
    """
    def __init__(self, master, count):
        tk.Frame.__init__(self, master, width=50, height=5, borderwidth=5, relief='groove')
        self.pack()
        tk.Label(self, text="Processing...").pack(padx=15, pady=10)

        self.progress = ttk.Progressbar(self, orient='horizontal', length=250, mode='determinate')
        self.progress.pack(padx=15, pady=10)
        self.progress['value'] = 0
        self.progress['maximum'] = count


class CustomToolbar(NavigationToolbar2TkAgg):
    """Creating a custom toolbar to be used with tkAgg
    """
    toolitems = filter(lambda x: x[0] != "Subplots", NavigationToolbar2TkAgg.toolitems)

    
if __name__ == '__main__':
    mw = QUIDDITMain("QUIDDIT version 2.0")
    mw.mainloop()
