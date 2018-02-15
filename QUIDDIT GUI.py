#-*- coding: utf-8 -*-
"""
Created on Thu Sep 07 12:31:42 2017

@author: ls13943
"""

import webbrowser
import tkinter as tk
import tkinter.filedialog as fd
from tkinter import ttk
import numpy as np

from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.cm as cm
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable


from tclwinbase_QUIDDITversion import *
import QUIDDIT_settings as settings
import QUIDDIT_utility as utility
import QUIDDIT_baseline
import QUIDDIT_main

QUIDDITversion = settings.version
STDBG = '#ececec'

all_toolitems = ('Home', 'Back', 'Forward', None,
                 'Pan', 'Zoom', 'Subplots', None,
                 'Save')


###############################################################################

class QUIDDITMain(TclWinBase):
    """Main window for QUIDDIT
    """

    def make_gui(self, title):
        """Die GUI im Haupt-Frame aufbauen
        """
        self.setwintitle(title)

        self.set_defaults()

        self.master.protocol('WM_DELETE_WINDOW', self.click_exit)
        self.master.bind('<Return>', self.click_ok)
        self.master.bind('<Escape>', self.click_exit)

        #Making the main menu
        menubar = tk.Menu(self)

        filemenuopt = {'Open...':self.file_open,
                       'Exit':self.click_exit}
        filemenu = self.make_menu(menubar, 'File', filemenuopt)
        filemenu.insert_separator(1)

        optmenuopt = {'settings':self.settings,
                      'Custom baseline':self.hello}
        self.make_menu(menubar, 'Settings', optmenuopt)

        baselinemenuopt = {'Correct baseline (default)':self.baseline,
                           'Correct with custom baseline':self.hello,
                           'Another option':self.hello}
        self.make_menu(menubar, 'Baseline', baselinemenuopt)

        procmenuopt = {'Process Data':self.process_data}
        self.make_menu(menubar, 'Process', procmenuopt)

        revmenuopt = {'Review fitting':self.review}
        self.make_menu(menubar, 'Review', revmenuopt)

        plotmenuopt = {'Plot single spectra':self.file_open,
                       'Plot line data':self.plot_ls,
                       'Plot map data':self.plot_map}
        self.make_menu(menubar, 'Plot', plotmenuopt)

        helpmenuopt = {'Read Manual':(lambda s=self: webbrowser.open(s.github_url)),
                       'About':self.about}
        helpmenu = self.make_menu(menubar, 'Help', helpmenuopt)
        helpmenu.insert_separator(1)

        self.master.config(menu=menubar)


##############################################################################

        row = 0

        self.fig = Figure(dpi=100)
        #self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.fig.gca().invert_xaxis()
        self.fig.suptitle('QUIDDIT')
        self.IIa_spec = np.loadtxt(settings.IIa_path, delimiter=',')
        self.ax.plot(self.IIa_spec[:, 0], self.IIa_spec[:, 1], 'k-')
        #img = mpimg.imread(r'C:/Users\ls13943\Dropbox\coding\QUIDDIT\QUIDDIT logo.gif')
        #self.ax.imshow(img)

        self.canvas = self.make_mplcanvas(self.fig, erow=row, ecol=0, cspan=4)

        self.message = tk.Text(self, state='disabled', relief=STDRELIEF)
        self.message.grid(row=row, column=4, columnspan=2, sticky=tk.NSEW, padx=5)
        scrl_bar = tk.Scrollbar(self, command=self.message.yview)
        scrl_bar.grid(row=row, column=6, sticky=tk.NSEW)
        self.message['yscrollcommand'] = scrl_bar.set

        row += 1

        toolbar_frame = tk.Frame(self)
        toolbar_frame.grid(row=row, column=1, sticky=tk.W)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, toolbar_frame)
        #self.toolbar = CustomToolbar(self.canvas, toolbar_frame)

        self.position = 0

        self.Back_button = self.makebutton(erow=row, ecol=0, caption='Previous',
                                           cmd=lambda s=self: s.plot_seq('prev'),
                                           state='disabled', padx=5, pady=5)
        self.Next_button = self.makebutton(erow=row, ecol=3, caption='Next',
                                           cmd=lambda s=self: s.plot_seq('next'),
                                           state='disabled', padx=5, pady=5, sticky='e')

        self.Histo_button = self.makebutton(erow=row, ecol=2, caption='Histogram',
                                           cmd=self.plot_histo,
                                           state='disabled', padx=5, pady=5, sticky='e')


        self.stats = self.maketext(lcol=4, lrow=row, erow=row, ecol=4,
                                   caption='stats: ', height=2, relief=STDRELIEF,
                                   state='disabled', background=STDBG)

        self.print_message(self.stats, self.file_count.get())

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
        self.logo.grid(padx=5, pady=5)

        about_msg = 'version {}\n(Laura Speich, 10/2017)\n\nlaura.speich@bristol.ac.uk'.format(str(QUIDDITversion))
        msg = tk.Label(self.toplevel, text=about_msg)
        msg.grid(padx=5, pady=5)
        Dismiss_button = tk.Button(self.toplevel, text='Dismiss',
                                   command=self.toplevel.destroy,
                                   height=1, width=6, default='active')
        Dismiss_button.grid(padx=5, pady=5)


    def baseline(self):
        """Baseline correct data
        """
        self.selected_files = fd.askopenfilenames(parent=self, initialdir=self.home,
                                                  title='Select CSV files',
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
                    self.update()
                    QUIDDIT_baseline.main(item, self.save_dir)
                loading.destroy()


    def clear_plot(self, fig):
        """Clear all axes in fig
        """
        for ax in fig.get_axes():
            fig.delaxes(ax)


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



    def file_open(self):
        """Open and plot spectrum
        """
        self.selected_items = fd.askopenfilenames(parent=self, initialdir=self.home,
                                                  title='Select spectra (CSV) files to display',
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
        self.set_entry_text(self.review_name, '[sample name] results')
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
        self.user_inp = (self.namevar.get(), self.resultvar.get(),
                         self.reviewvar.get(), self.agevar.get())

        return self.user_inp


    def loaded(self):
        """Das Window wurde aufgebaut
        """
        self.print_message(self.message, 'Welcome to QUIDDIT ver. 2.0\n')
        self.std = settings.std


    def on_input(self, event):
        """Get sample name and use it to suggest file names
        """
        inp = self.namevar.get()
        self.result_name.delete('0', 'end')
        self.result_name.insert('0', inp+' results')
        self.review_name.delete('0', 'end')
        self.review_name.insert('0', inp+' review')

    def plot_histo(self):
        self.toplevel = QUIDDITToplevel('Histogram')
        
        row = 0
        
        self.fig = Figure(dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.fig.suptitle(self.map_title)
        self.ax.hist(self.plot_item[~np.isnan(self.plot_item)], bins=100)
        
        self.canvas = self.toplevel.make_mplcanvas(self.fig, erow=row, ecol=0,
                                                   rspan=2)

        self.min = self.toplevel.make_double_entry(lcol=1, lrow=row, ecol=2, erow=row,
                                          caption='min')

        row += 1

        self.max = self.toplevel.make_double_entry(lcol=1, lrow=row, ecol=2, erow=row,
                                          caption='max')

        row += 1

        self.toplevel.makebutton(erow=row, ecol=1, cspan=2,
                                 padx=5, pady=5, sticky='nsew',
                                 caption='Redo map', cmd=self.hello)

        toolbar_frame = tk.Frame(self.toplevel)
        toolbar_frame.grid(row=row, column=0, sticky=tk.W, columnspan=3)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, toolbar_frame)
    

    def plot_ls(self):
        """Method to plot results of linescans.
        """
        self.res_file = fd.askopenfilename(parent=self, initialdir=self.home,
                                           title='Select results file',
                                           filetypes=(('CSV', '*.CSV'), ('CSV', '*.csv')))
        self.clear_plot(self.fig)
        #self.fig.set_size_inches(8, 7)
        #self.canvas.config(width=8, height=7)
        self.fig.suptitle('')

        res = np.loadtxt(self.res_file, delimiter=',', dtype=utility.results_dtype, skiprows=2)
        points = []
        for item in res['name']:
            points.append(float(item.split()[-1][-8:-4]))

        self.ax0 = self.fig.add_subplot(3, 3, 1)
        self.ax0.text(0, 0.5, self.res_file.split('/')[-1], bbox={'facecolor':'white', 'pad':10})
        self.ax0.axis('off')

        self.ax1 = self.fig.add_subplot(3, 3, 2)
        self.ax1.plot(res['[NB]']/79.4, res['p_area_ana'], 'k.')
        self.ax1.plot((0, 100), (0, 6400), 'k-')
        self.ax1.set(xlim=(0, 10), ylim=(0, 640),
                     xlabel='$\mathregular{\mu_B}$', ylabel="I(B') $\mathregular{(cm^{-2})}$")

        self.ax2 = self.fig.add_subplot(3, 3, 3)
        self.ax2.plot(points, res['H_area_ana'], 'k.')
        self.ax2.set(ylim=(0, None),
                     ylabel='I(3107) $\mathregular{(cm^{-2})}$',
                     xticklabels=[])

        self.ax3 = self.fig.add_subplot(3, 3, 4)
        self.ax3.plot(points, res['[NA]'], '.', label='$\mathregular{[N_A]}$')
        self.ax3.plot(points, res['[NB]'], '.', label='$\mathregular{[N_B]}$')
        self.ax3.plot(points, res['[NT]'], '.', label='$\mathregular{[N_t]}$')
        self.ax3.set(ylim=(0, None),
                     ylabel='concentration (ppm)',
                     xticklabels=[])
        self.ax3.legend(loc='best')

        self.ax4 = self.fig.add_subplot(3, 3, 5)
        self.ax4.plot(points, (res['[NB]']/(res['[NA]']+res['[NB]'])), 'k.')
        self.ax4.set(ylim=(0, 1),
                     ylabel='$\mathregular{[N_B]/[N_T]}$',
                     xticklabels=[])

        self.ax5 = self.fig.add_subplot(3, 3, 6)
        self.ax5.plot(points, res['T'], 'k.')
        self.ax5.set(ylim=(1050, 1350),
                     ylabel='$\mathregular{T_N (^{\circ}C)}$',
                     xticklabels=[])

        self.ax6 = self.fig.add_subplot(3, 3, 7)
        self.ax6.plot(points, res['p_area_ana'], 'k.')
        self.ax6.set(ylim=(0, None),
                     xlabel='#', ylabel="I(B') $\mathregular{(cm^{-2})}$")

        self.ax7 = self.fig.add_subplot(3, 3, 8)
        P0 = (res['[NB]']/79.4)*64
        deg = (1-(res['p_area_ana']/P0))
        self.ax7.plot(points, deg, 'k.')
        self.ax7.set(ylim=(0, 1),
                     xlabel='#', ylabel="$\mathregular{1-(I(B')/I(B')_0)}$")

        self.ax8 = self.fig.add_subplot(3, 3, 9)
        lnstuff = np.log(np.log(P0/res['p_area_ana'])/(1000*1e6*365*24*60*60))
        TP = ((88446)/(19.687-lnstuff)) -273
        self.ax8.plot(points, TP, 'k.')
        self.ax8.set(ylim=(1050, 1350),
                     xlabel='#', ylabel='$\mathregular{T_P (^{\circ}C)}$')

        self.fig.tight_layout(pad=0.7, w_pad=0, h_pad=0)
        self.canvas.draw()


    def plot_map(self):
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


    def plot_seq(self, seq):
        """plots spectra within a sequence.
        plotmode=='single': plot single spectrum
        plotmode=='review': plot spectrum and review elements
        plotmode=='map': plot maps
        """
        self.clear_plot(self.fig)

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
            self.ax = self.fig.add_subplot(1, 1, 1)
            self.plot_spec((self.selected_items[self.index]),
                           self.fig, self.ax)

        elif self.plotmode.get() == 'review':
            self.fig.suptitle('Review')
            for i in range(3):
                self.ax = self.fig.add_subplot(3, 1, i+1)
                self.plot_spec(self.selected_items[self.index],
                               self.fig, self.ax, xlabel='')
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
                p_peak_area = np.round(utility.peak_area(*p_params[:5]), 1)

                H_spec = utility.spectrum_slice(spec, 3000, 3200)
                H_wav = H_spec[:, 0]
                H_params = (rev['H_pos'], rev['H_I'], rev['H_HWHM_l'], rev['H_HWHM_r'], rev['H_sigma'])
                H_bg_params = (rev['H_bg_a'], rev['H_bg_b'], rev['H_bg_c'], rev['H_bg_d'])
                H_fit = utility.pseudovoigt_fit(H_wav, *H_params) + np.polyval(H_bg_params, H_wav)
                H_peak_area = np.round(utility.peak_area(*H_params), 1)

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
                                       '3107 peak area: {0:.2f} cm-2'.format(H_peak_area))

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
            self.Histo_button['state'] = 'normal'
            x = []
            y = []
            res = self.read_results(self.res_file)
            for item in res['name']:
                xy = item.decode().split('/')[-1]
                x.append(float(xy.split(' ')[0][2:]))
                y.append(float(xy.split(' ')[1][1:-5]))

            extent = (min(x), max(x), min(y), max(y))
            resolution = 2000j
            cols = cm.jet

            grid_x, grid_y = np.mgrid[extent[0]:extent[1]:resolution,
                                      extent[2]:extent[3]:resolution]

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

            self.plot_item = maps[plots[self.index]]
            self.map_title = plots[self.index]

            if plots[self.index] == '$[N_B]/[N_T]$':
                self.ax.set()

            map_grid = utility.make_2dgrid(x, y, grid_x, grid_y, self.plot_item)
            self.fig.suptitle(self.map_title)
            self.ax = self.fig.add_subplot(1, 1, 1)
            img = self.ax.imshow(map_grid, origin='lower', extent=extent,
                                 cmap=cols, clim=clims[plots[self.index]])
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.3)
            self.cbar = self.fig.colorbar(img, cax=cax)

        else:
            self.print_message(self.message, 'Error: Unknown plotmode.')

        self.canvas.draw()


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
                                                  title='Select CSV files',
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

                self.results[i] = curr_res
                self.review[i] = curr_rev

                loading.progress['value'] += 1
                #loading.update()
                self.update()
                i += 1

            loading.destroy()

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
        #self.file_count = tk.StringVar(value='')
        self.file_count = self.getvar('')
        self.github_url = 'https://github.com/LauraSp/QUIDDIT'
        self.namevar = self.getvar('')
        self.resultvar = self.getvar('')
        self.reviewvar = self.getvar('')
        self.agevar = self.getvar(2900)
        #self.selected_files = self.getvar('')
        #self.selected_items = self.getvar('')
        self.plotmode = self.getvar('')


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


class QUIDDITToplevel(tk.Toplevel, TclWinBase):
    """Making a toplevel window that inherits from TclWinBase
    """
    def __init__(self, title):
        super().__init__()
        self.toplevel = tk.Toplevel()
        #self.toplevel.title(title)
        self.protocol("WM_DELETE_WINDOW", self.toplevel.destroy)


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