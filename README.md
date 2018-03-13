# QUIDDIT
Quantification of Infrared active Defects in Diamond and Inferred Temperatures

license: 

## Download
All necessary files can be downloaded by visiting the following website:
https://github.com/LauraSp/QUIDDIT

A full manual is available at this address as well.

## Installation instructions
### Install Python
In order to run QUIDDIT, you will need a working version of Python 3.

I recommend installing an integrated development environment (IDE) that includes the most commonly used libraries. To run QUIDDIT, you will need:
* SciPy (https://www.scipy.org/)
* NumPy (http://www.numpy.org/)
* matplotlib (https://matplotlib.org/)
* Tkinter (https://docs.python.org/2/library/tk.html)
* webbrowser (https://docs.python.org/2/library/webbrowser.html)

All of these are part of the most common IDE for Python. The instructions in the manual were created for use with Spyder (which is part of the Anaconda package), so I recommend installing Anaconda and running scripts with Spyder for users not familiar with coding.

To install Anaconda, visit https://www.anaconda.com/download/

Chose your operating system and follow the instructions.

### Install QUIDDIT
Download the QUIDDIT package from https://github.com/LauraSp/QUIDDIT, unpack if necessary.
Create a directory on your computer to store spectral data and processing results in, such as
C:\FTIR

### Running QUIDDIT
To run QUIDDIT, open your python IDE of choice and find and run the "QUIDDIT GUI" file in the downloaded repository (you may need to unzip the files first). The GUI can also be run using IDLE (a very basic standard Python IDE) or from the command line.

## Known Bugs and Issues
This section provides an overview over known issues with QUIDDIT in order of priority. The author is working on resolving them but no guarantee can be given at what point they will be fixed.

* Issues in fitting the platelet peak area (fit overestimates height of the 1405 cm-1 peak)
* Results of spectra that contain the C-centre currently can't be plotted
* When plotting map data, "Redo map" button in "Histogram" window only works the second time it is clicked
* Empty, additional pop-up window prompted with "Settings" and "User Input" (in data processing)

## Contact:
For further questions or suggestions, please contact

Laura Speich: laura.speich@bristol.ac.uk
