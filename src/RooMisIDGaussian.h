/*****************************************************************************
 * Project: pdf for Gaussian with mis ID invariant mass correction           *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 01-2016                                                            *
 *****************************************************************************/

#ifndef ROOMISIDGAUSSIAN
#define ROOMISIDGAUSSIAN

#include "RooAbsPdf.h"
#include "RooAbsProxy.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooMomentumFractionPdf.h"
 
class RooMisIDGaussian : public RooAbsPdf {
public:
  RooMisIDGaussian() {}; 
  
  // two dimensional (mass and momentum fraction)
  RooMisIDGaussian(const char *name, const char *title,
          RooAbsReal& _m,
          RooAbsReal& _mean,
          RooAbsReal& _sigma,
          RooAbsReal& _dm2,
          RooMomentumFractionPdf& _xpdf);
  
  RooMisIDGaussian(const char *name, const char *title,
          RooAbsReal& _m,
          RooAbsReal& _mean,
          RooAbsReal& _sigma,
          RooAbsReal& _dm2,
          RooAbsReal& _x,
          RooAbsPdf& _xpdf);

  RooMisIDGaussian(const char *name, const char *title,
          RooAbsReal& _m,
          RooAbsReal& _mean,
          RooAbsReal& _sigma,
          RooAbsReal& _dm2,
          RooAbsReal& _x);
  
  RooMisIDGaussian(const char *name, const char *title,
          RooAbsReal& _m,
          RooAbsReal& _mean,
          RooAbsReal& _sigma,
          RooAbsReal& _dm2,
	      RooAbsReal& _x,
          RooAbsReal& _xmin,
          RooAbsReal& _dx,
          RooArgList& _slopelist,
	      RooArgList& _coefslist);

  // one dimensional (mass only)
  RooMisIDGaussian(const char *name, const char *title,
          RooAbsReal& _m,
          RooAbsReal& _mean,
          RooAbsReal& _sigma,
          RooAbsReal& _dm2,
          RooAbsReal& _power,
          RooAbsReal& _xmin,
          RooAbsReal& _xmax,
          RooAbsReal& _dx);

  RooMisIDGaussian(const char *name, const char *title,
          RooAbsReal& _m,
          RooAbsReal& _mean,
          RooAbsReal& _sigma,
          RooAbsReal& _dm2,
          RooAbsReal& _xmin,
          RooAbsReal& _dx,
          RooArgList& _slopelist,
	      RooArgList& _coefslist,
          RooAbsReal& _xmax);

  RooMisIDGaussian(const RooMisIDGaussian& other, const char* name=0) ;
  
  virtual TObject* clone(const char* newname) const { return new RooMisIDGaussian(*this,newname); }
  inline virtual ~RooMisIDGaussian() { }

  // analytical integral functions for RooFit
  Int_t	getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

  // function to set integrals to numerical
  inline void setNumericalIntegral(const Bool_t flag_nm = true) { flag_num_int = flag_nm; } ;

protected:
  RooRealProxy m ;
  RooRealProxy mean ;
  RooRealProxy sigma ;
  RooRealProxy dm2 ;
  RooRealProxy x ;
  RooRealProxy xmin ;
  RooRealProxy xmax ;
  RooRealProxy dx ;
  RooRealProxy power ;
  RooRealProxy abspdf ; 
  RooListProxy slopelist ;
  RooListProxy coefslist ;

  Int_t slopeSize ;
  TIterator* _slopeIter ;
  TIterator* _coefsIter ;

  // flag for conditional x 
  Bool_t flag_conditional_x ;

  // flag for multidimensional (m,x) dependence
  Bool_t flag_2D ;

  // flag for numerical integrals
  Bool_t flag_num_int ;

  // flag for use of constructor with abstract pdf 
  Bool_t flag_abspdf ;

  // Error function approximator
  static Double_t ApproxErf(const double arg);

  // primitive functions 
  Double_t primitiveExpoGauss2D(const double mval, 
        const double xval, const double meanval, const double sigval, 
        const double dm2val, const double cval, const double aval) const;
  Double_t primitiveExpoGauss1D(const double mval, 
        const double xval, const double meanval, const double sigval, 
        const double dm2val, const double cval, const double aval) const;
  Double_t primitiveSplineGauss2D(const double mval, 
        const double xval, const double meanval, const double sigval, 
        const double dm2val, const double dxval, const double xminval) const;
  Double_t primitiveSplineGauss1D(const double mval, 
        const double xval, const double meanval, const double sigval, 
        const double dm2val, const double dxval, const double xminval) const; 

  // evaluation helper functions
  Double_t getValExponentials(const double xval) const;
  Double_t getValSpline(const double xval) const;
  Double_t getValGaussXintegrated(const double mval) const;

  // integration over x
  Double_t getIntExponentialsGaussOverX(const double x0, const double x1, const double mval) const;
  Double_t getIntSplineGaussOverX(const double x0, const double x1, const double mval) const;

  // integration over m
  Double_t getIntGaussXintegrated(const double m0, const double m1) const;
  Double_t getIntGaussNotXintegrated(const double m0, const double m1, const double xval) const;

  // multi dimensional integrals
  Double_t getMultiIntExponentialsGauss(const double x0, const double x1,
          const double m0, const double m1) const;
  Double_t getMultiIntSplineGauss(const double x0, const double x1,
          const double m0, const double m1) const;

  // general evaluation function
  Double_t evaluate() const ;

private:
  ClassDef(RooMisIDGaussian,1) 
};
 
#endif
