This programme can automatically process the files from the Tektronix 2024b oscilloscope (providing the raw .csv file has had the setup removed, and the measurements have been moved to the first two columns). 

The data file will be name "TEKnnnn.CSV". The files for the preamp data should be an even number, and the square wave should be odd 

The preamp data should be taken with a square wave of a function generator that charges up a 3.3 pF capacitor. 

The oscilloscope needs to be triggered on the falling edge of the square wave, and should be zoomed in so that the time division is 1us. This is to ensure that the preamp peak approximates a square waveform. 

The programme should be run as:

"./preamp-fitting command1 command2 command3"

Command1 - Can either be "-a", or "-f". -a will analyse the raw data to append the signal charge and error, as well as the output voltage and error to a file name "test.csv". Ensure that the file is either empty if starting a new analysis, or if it is the correct set of data already in it if continuing an analysis. -f will fit the aforemention "test.csv" file, and output the gradient and its error, the intercept and its error, and the chi2. It will also output a pdf of the fit and calibration curve.

Command2 - Will be the file name of either the preamp data (including the foldername) if command1 is -a, or the name of the file containing the calibration curve if command 1 is -f

command3 - Will be the name of the square wave data if command1 is -a. Will not be used if command1 is -f
