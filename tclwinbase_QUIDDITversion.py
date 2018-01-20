"""My Version of using Tkinter
"""
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

STDRELIEF = tk.FLAT

class TclWinBaseUsageException(BaseException):
    """Exception für Nutzungsfehler dieses Moduls
    """
    def __init__(self, arg):
        self.args = arg


class TclWinBase(tk.Frame):
    """base class for tcl driven windows
    inherit this for your own windows
    """

    def mainloop(self):
        """ Enter the main loop"""
        self.root.mainloop()

    def __init__(self, title):
        #NoDefaultRoot()
        self.root = tk.Tk() 
        super().__init__(self.root, padx=5, pady=5)
        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_rowconfigure(0, weight=1)
        self.root.resizable(True, True)
        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)

        self.make_gui(title)
        self.loaded()

    def setwintitle(self, title):
        """Sets the title of your main window
        """
        self.root.title(title)

    def loaded(self):
        raise TclWinBaseUsageException("Override me! Always override loaded method")

    def make_gui(self, title):
        raise TclWinBaseUsageException("Override me! Always override make_gui method")
        
    def make_menu(self, menubar, title, itemlib):     
        submenu = tk.Menu(menubar, tearoff=0)
        for item in itemlib:
            submenu.add_command(label=item, command=itemlib[item])
            
        menubar.add_cascade(label=title, menu=submenu)
        
        return submenu
       

    def maketext(self, lcol=0, lrow=0, erow=0, ecol=1, cspan=1, rspan=1, caption='', width=None, **options):
        """create a multiple single line text widget with a label/caption in another column
        """
        tk.Label(self, text=caption).grid(row=lrow, column=lcol, columnspan=cspan, rowspan=rspan, sticky=tk.N + tk.E)
        entry = tk.Text(self, **options)
        if width:
            entry.config(width=width)
    
        entry.grid(row=erow, column=ecol, sticky=tk.W)
        return entry

    def makeentry(self, lcol=0, lrow=0, erow=0, ecol=1, caption='', width=None, **options):
        """create a single line text entry widget with a label"""
        tk.Label(self, text=caption).grid(row=lrow, column=lcol, sticky=tk.E)
        entry = tk.Entry(self, relief=STDRELIEF, **options)
        if width:
            entry.config(width=width)
    
        entry.grid(row=erow, column=ecol, sticky=tk.W)
        return entry

    def make_double_entry(self, lcol=0, lrow=0, erow=0, ecol=1, caption='', width=None, **options):
        """create a single line text for number entry widget with a label
        """
        tk.Label(self, text=caption).grid(row=lrow, column=lcol, sticky=tk.E)
        entry = ValidateDoubleEntry(self, relief=STDRELIEF, **options)
        if width:
            entry.config(width=width)
    
        entry.grid(row=erow, column=ecol, sticky=tk.W)
        return entry

    def make_int_entry(self, lcol=0, lrow=0, erow=0, ecol=1, caption='', width=None, **options):
        """create a single line text for number entry widget with a label
        """
        tk.Label(self, text=caption).grid(row=lrow, column=lcol, sticky=tk.E)
        entry = ValidateIntegerEntry(self, relief=STDRELIEF, **options)
        if width:
            entry.config(width=width)
    
        entry.grid(row=erow, column=ecol, sticky=tk.W)
        return entry

    def set_entry_text(self, entry, text):
        """Set text in text entry to a given text"""
        entry.delete(0, tk.END)
        entry.insert(tk.END, text)

    def makecheck(self, ecol=0, erow=0, caption='', **options):
        """create a checkbox with a label"""
        cb = tk.Checkbutton(self, text=caption, **options)
        cb.grid(row=erow, column=ecol, sticky=tk.W)
        return cb

    def makebutton(self, erow=0, ecol=0, cspan=1, rspan=1, caption='Button', 
                   width=None, cmd=None, sticky=tk.W, **options):
        """create a button widget"""
        bu = tk.Button(self,
                        text=caption,
                        width=width,
                        command=cmd, 
                        **options)
        
        bu.grid(row=erow, column=ecol, columnspan=cspan, rowspan=rspan, sticky=sticky)

        return bu

    def makecanvas(self, erow=0, ecol=0, rspan=1, cspan=1, sticky=tk.NSEW, **options):
        """create a canvas widget"""
        ca = tk.Canvas(self, **options)
        ca.grid(row=erow, column=ecol,
                columnspan=cspan, rowspan=rspan, sticky=sticky)

        return ca
    
    def make_mplcanvas(self, fig, erow=0, ecol=0, rspan=1, cspan=1, sticky=tk.NSEW, **options):
        canvas = FigureCanvasTkAgg(fig, master=self, **options)
        mpl_canvas = canvas.get_tk_widget()
        mpl_canvas.grid(row=erow, column=ecol, columnspan=cspan, rowspan=rspan, sticky=sticky)
        
        return canvas

    def makelist(self, lcol=0, lrow=0, erow=0, ecol=1, caption='', width=None,
                 scrollvert=True, scrollhor=False,
                 **options):
        """create a list widget in the current window
        """
        
        tk.Label(self, text=caption).grid(row=lrow, column=lcol, sticky=tk.N + tk.E)
        
        if scrollvert == True:
            yScroll = tk.Scrollbar(self, orient=tk.VERTICAL)
            yScroll.grid(row=erow, column=ecol+1, sticky=tk.N+tk.S)
            
        if scrollhor == True:
            xScroll = tk.Scrollbar(self, orient=tk.HORIZONTAL)
            xScroll.grid(row=erow+1, column=ecol, sticky=tk.E+tk.W)

        lst = tk.Listbox(self, **options)
        lst.grid(row=erow, column=ecol)
        
        if scrollvert == True:
            lst.config(yscrollcommand=yScroll.set)
            yScroll['command'] = lst.yview

        if scrollhor == True:
            lst.config(xscrollcommand=xScroll.set)
            xScroll['command'] = lst.xview

        if width:
            lst.config(width=width)

        return lst

    def getvar(self, defval):
        """gets a new tkinter variable to be used for binding to entry widgets
        """
        t = type(defval)

        if t == str:
            answ = tk.StringVar()
            answ.set(defval)
        elif t == int:
            answ = tk.IntVar()
            answ.set(defval)
        elif t == float:
            answ = tk.DoubleVar()
            answ.set(defval)
        else:
            answ = tk.StringVar()

        return answ


class ValidateDoubleEntry():
    """Special entry widget for editing of double numbers
    """
    def __init__(self, parent, **options):
        validate_number_cmd = parent.register(self.validate_number)
        self.entry = tk.Entry(parent,
                               validate='all',
                               validatecommand=(validate_number_cmd, '%d', '%i', '%S'),
                               **options)
    
    def config(self, **options):
        self.entry.config(**options)

    def grid(self, **options):
        self.entry.grid(**options)

    def validate_number(self, d, i, s):
        if s == '':
            return True

        if s.isdigit() or s=='.' or (i=='0' and s=='-'):
            return True

        return False


class ValidateIntegerEntry():
    """special entry type widget for editing integer values
    """
    def __init__(self, parent, **options):
        validate_number_cmd = parent.register(self.validate_number)
        self.entry = tk.Entry(parent,
                               validate='all',
                               validatecommand=(validate_number_cmd, '%d', '%i', '%S'),
                               **options)
    
    def config(self, **options):
        self.entry.config(**options)

    def grid(self, **options):
        self.entry.grid(**options)

    def validate_number(self, d, i, s):
        if s == '':
            return True

        if s.isdigit() or (i=='0' and s=='-'):
            return True

        return False

