# RooMisID
A RooFit package for fitting misidentified backgrounds

###Contains the following classes
  * __RooMisIDGaussian__  
    Main class for a probability density function (PDF) for invariant mass for a
    decay with one particle misidentified, based on a gaussian mass resolution where the mean depedends on the kinematical variable called "momentum fraction"
  * __RooMomentumFractionPdf__  
    Class for PDF for the momentum fraction, in its core a sum of exponentials with a 3rd order spline (polynomial) at its start (low x). 
  * __RooMisIDBkg__
    Class for PDF for invariant mass for a decay with one particle misidentified, without resolution.
    Suggestion is to use it with your own resolution pdf and apply a FFT convolution (RooFFTConvPdf).
  * __RooInverseGaussian__  
    Class for an inverse Gaussian PDF
  * __RooPowerLawApproxVar__  
    Class used in RooMisIDGaussian and RooMomentumFractionPDF to approximate a power law with exponentials

###Build instructions
  1. make a folder to build in (recommended):  
   ```
   cd <directory_where_package_lives>/RooMisID
   ```  
   ```
   mkdir build
   ```  
   ```
   cd build  
   ```
  2. use cmake to set up environment and build here (ensure you have cmake!):  
   ```
   cmake ../src
   ```  
   ```
   make  
   ```
  3. the library lives in:  
   ```
   <directory_where_package_lives>/RooMisID/build/libRooMisID.so  
   ```
