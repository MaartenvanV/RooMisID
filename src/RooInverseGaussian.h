/*****************************************************************************
 * Project: inverse gaussian pdf                                             *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 02-2016                                                            *
 *****************************************************************************/

#ifndef ROOINVERSEGAUSSIAN    
#define ROOINVERSEGAUSSIAN    

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"

class RooInverseGaussian : public RooAbsPdf {
    public:
        RooInverseGaussian() {}; 

        RooInverseGaussian(const char *name, const char *title,
                RooAbsReal& _x,
                RooAbsReal& _xmin,
                RooAbsReal& _lambda,
                RooAbsReal& _mu);

        RooInverseGaussian(const RooInverseGaussian& other, const char* name=0) ;

        virtual TObject* clone(const char* newname) const { return new RooInverseGaussian(*this,newname); }
        inline virtual ~RooInverseGaussian() { }

        // integration functions
        Int_t	getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
        Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;
        
        // generator functions
        Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
        void generateEvent(Int_t code);

    protected:
        RooRealProxy m_x ;
        RooRealProxy m_xmin ;
        RooRealProxy m_lambda ;
        RooRealProxy m_mu ;

        Double_t evaluate() const ;

    private:
        ClassDef(RooInverseGaussian,1) 
};

#endif
