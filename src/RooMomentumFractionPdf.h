/*****************************************************************************
 * Project: Momentum fraction pdf for mis ID invariant mass correction       *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 01-2016                                                            *
 *****************************************************************************/

#ifndef ROOMOMENTUMFRACTIONPDF
#define ROOMOMENTUMFRACTIONPDF

#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
 
class RooMomentumFractionPdf : public RooAbsPdf {
public:
  RooMomentumFractionPdf() {}; 
 
  RooMomentumFractionPdf(const char *name, const char *title,
          RooAbsReal& _x,
          RooAbsReal& _xmin,
          RooAbsReal& _dx,
          RooAbsReal& _power);
  
  RooMomentumFractionPdf(const char *name, const char *title,
	      RooAbsReal& _x,
          RooAbsReal& _xmin,
          RooAbsReal& _dx,
          RooArgList& _slopelist,
	      RooArgList& _coefslist);

  RooMomentumFractionPdf(const RooMomentumFractionPdf& other, const char* name=0) ;
  
  virtual TObject* clone(const char* newname) const { return new RooMomentumFractionPdf(*this,newname); }
  inline virtual ~RooMomentumFractionPdf() { }

  // integration functions
  Int_t	getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

  // parameters retrieval function
  RooArgList* parseDependents() const;

protected:
  RooRealProxy m_x ;
  RooRealProxy m_xmin ;
  RooRealProxy m_dx ;
  RooRealProxy m_power ;
  RooListProxy m_slopelist ;
  RooListProxy m_coefslist ;

  Int_t m_slopeSize ;
  TIterator* m_slopeIter ;
  TIterator* m_coefsIter ;

  // helper functions
  Double_t getValExponentials(const double xval) const;
  Double_t getValSpline(const double xval) const;
  Double_t getIntExponentials(const double x0, const double x1) const;
  Double_t getIntSpline(const double x0, const double x1) const;

  Double_t evaluate() const ;

private:
  ClassDef(RooMomentumFractionPdf,1) 
};
 
#endif
