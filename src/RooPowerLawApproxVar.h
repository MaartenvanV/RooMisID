/*****************************************************************************
 * Project: Helper class for coefficients and slopes of exponentials         *
 *          approximating a power law                                        *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 01-2016                                                            *
 *****************************************************************************/

#ifndef ROOPOWERLAWAPPROXVAR
#define ROOPOWERLAWAPPROXVAR

#include "RooRealProxy.h"
#include "RooAbsReal.h"
 
class RooPowerLawApproxVar : public RooAbsReal {
public:
  // flags
  typedef enum {
      Coeff = 1,
      Slope = 2,
  } Flags;

  // constructors
  RooPowerLawApproxVar() {};

  RooPowerLawApproxVar(const char* name, const char *title,
          RooAbsReal& _power,
          const Int_t _index,
          const Double_t _beta_approx,
          const Int_t _n_expos,
          Flags _flags);

  RooPowerLawApproxVar(const RooPowerLawApproxVar& other, const char* name=0) ;

  virtual TObject* clone(const char* newname) const { return new RooPowerLawApproxVar(*this,newname); }
  inline virtual ~RooPowerLawApproxVar() { }

protected:
  RooRealProxy m_power ;
  Double_t m_beta_approx ;
  Int_t m_n_expos ;
  Int_t m_index ;
  Flags m_flags ;

  Double_t evaluate() const ;

private:
  ClassDef(RooPowerLawApproxVar,1) 
};

#endif
