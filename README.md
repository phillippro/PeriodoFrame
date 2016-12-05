# PeriodoFrame
A framework of periodograms to disentangle periodic signals from correlated noise

Since these codes are wirtten in R, R and relevant packages are required to run the code. The "minpack.lm" R package is required to calculate the BFP. The "fields" and "magicaxis" packages are optional, and are needed to make figures.  

The test data is put in the data/ directory. And the results are put in the results/ directory. One can change these directories. In the data/ directory, "...calibration.dat" files are calibration data, "..._TERRA_HARPS.dat" are TERRA-reduced HARPS data including JD-2400000, RV, eRV, BIS, FWHM and S-index from left to right columns. "..._TERRA_XAP_dRVs.dat" are differential RVs derived from XAP aperture data sets. In each file, the first column is the time, and the other columns are XAPj-i where j={2,3,...,X} and i=j-1.  

The functions needed to make PeriodoFrame are put in periodoframe.R while the functions used to make other periodograms are put into periodogrames.R . In specific, the MLP(...) and BFP(...) functions in periodoframe.R are used to create MLP and BFP respectively. The MP(...) function is used to make the moving periodogram. 

The other .R files are just examples for the usage of PeriodoFrame. The "make_periodoframe.R" file is used to make the BFP and MLP. The "MP.R" file is to make the moving periodogram. Both files call "prepare_data.R" which is to set parameters and to load the data. The reader is refered to the upcoming paper titled "PeriodoFrame: Disentangling periodic singals from correlated noise in a periodogram framework" for more details.  
