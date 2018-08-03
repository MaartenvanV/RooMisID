/*****************************************************************************
 * Project: pdf for mis ID invariant mass correction, without resolution     *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 08-2018                                                            *
 *****************************************************************************/

#ifndef ROOMISIDBKG
#define ROOMISIDBKG

#include "RooAbsPdf.h"
#include "RooAbsProxy.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
 
class RooMisIDBkg : public RooAbsPdf {
public:
  RooMisIDBkg() {}; 
 
  RooMisIDBkg(const char *name, const char *title,
          RooAbsReal& _m,
          RooAbsReal& _mean,
          RooAbsReal& _dm2,
          RooAbsReal& _power,
          RooAbsReal& _xmin,
          RooAbsReal& _xmax,
          RooAbsReal& _dx);

  RooMisIDBkg(const char *name, const char *title,
          RooAbsReal& _m,
          RooAbsReal& _mean,
          RooAbsReal& _dm2,
          RooAbsReal& _xmin,
          RooAbsReal& _xmax,
          RooAbsReal& _dx,
          RooArgList& _slopelist,
	      RooArgList& _coefslist);

  RooMisIDBkg(const RooMisIDBkg& other, const char* name=0) ;
  
  virtual TObject* clone(const char* newname) const { return new RooMisIDBkg(*this,newname); }
  inline virtual ~RooMisIDBkg() { }

  // integration functions
  Int_t	getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

protected:
  RooRealProxy m ;
  RooRealProxy mean ;
  RooRealProxy dm2 ;
  RooRealProxy x ;
  RooRealProxy xmin ;
  RooRealProxy xmax ;
  RooRealProxy dx ;
  RooRealProxy power ;
  RooListProxy slopelist ;
  RooListProxy coefslist ;

  Int_t slopeSize ;
  TIterator* _slopeIter ;
  TIterator* _coefsIter ;

  // flag for use of constructor with abstract pdf 
  Bool_t flag_abspdf ;

  // helper functions
  Double_t getValExponentials(const double xval) const;
  Double_t getValSpline(const double xval) const;

  // general evaluation function
  Double_t evaluate() const ;

private:
  ClassDef(RooMisIDBkg,1) 
};
 
#endif
