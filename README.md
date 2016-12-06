# PeriodoFrame
A framework of periodograms to disentangle periodic signals from correlated noise

1.Installation and usage
The user just download the code and put it in a directory. Since these codes are wirtten in R, R and relevant packages are required to run the code. The "minpack.lm" R package is required to calculate the BFP. The "fields" and "magicaxis" packages are optional, and are needed to make figures. Then the user can write their own code refereing to the examples provided by "make_periodoframe.R", "MP.R" and "prepare_data.R". The author can also use the functions directly to make BFP and MLP. 
For example, the BFP can be made using the following function:

bf <- BFP(t,y,dy,Nma=Nma,NI=NI,Indices=Indices,ofac=ofac,opt.type=opt.type,model.type=model.type,fmax=fmax,tol=tol)

where t, y and dy are the times, y values and the corresponding errors. Nma and NI are the numbers of moving average components and noise indices which are given by "Indices". ofac is the over sampling factor, and opt.type can be "nl" or "sl", and is the method used to optimize parameters. "sl" is suggested because it is more efficient and is as reliable as "nl". "model.type" is a parameter to specify whether to compare models. If model.type='auto', model comparison is done using the BFP to select the optimal noise model. If model.type='man', BFP will use Nma and NI manually provided by the user. "fmax" is the maximum frequency for frequency sampling, and "tol" is the tolorence or prevision required to optimize parameters (or maximize the likelihood) numerically. 

The main output of the BFP function is the logarithmic Bayes factor (logBF) for a sample of periods (P). Then the user can make a simple BFP using plot(P,logBF,log='x')

The usage of the functions of MLP and MP in periodoframe.R is similar, the user is refered to "make_periodoframe.R" and "MP.R" for examples. 


2.Test data
The test data is put in the data/ directory. And the results are put in the results/ directory. One can change these directories. In the data/ directory, "...calibration.dat" files are calibration data, "..._TERRA_HARPS.dat" are TERRA-reduced HARPS data including JD-2400000, RV, eRV, BIS, FWHM and S-index from left to right columns. "..._TERRA_XAP_dRVs.dat" are differential RVs derived from XAP aperture data sets. In each file, the first column is the time, and the other columns are XAPj-i where j={2,3,...,X} and i=j-1.  


3.Content of code files
The functions needed to make PeriodoFrame are put in periodoframe.R while the functions used to make other periodograms are put into periodogrames.R . In specific, the MLP(...) and BFP(...) functions in periodoframe.R are used to create MLP and BFP respectively. The MP(...) function is used to make the moving periodogram. 

The other .R files are just examples for the usage of PeriodoFrame. The "make_periodoframe.R" file is used to make the BFP and MLP. The "MP.R" file is to make the moving periodogram. Both files call "prepare_data.R" which is to set parameters and to load the data. The reader is refered to the upcoming paper titled "PeriodoFrame: Disentangling periodic singals from correlated noise in a periodogram framework" for more details.  
