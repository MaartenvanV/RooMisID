/*****************************************************************************
 * Project: inverse gaussian pdf                                             *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 02-2016                                                            *
 *****************************************************************************/

#include "RooFit.h"
#include "Riostream.h" 

#include <cassert>
#include <cmath>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "RooRandom.h"

#include "RooInverseGaussian.h" 

static const Double_t twopi = 6.283185307179586;
static const Double_t roottwo = 1.4142135623730951;

// constructors
RooInverseGaussian::RooInverseGaussian(const char *name, const char *title, 
            RooAbsReal& _x,
            RooAbsReal& _xmin,
            RooAbsReal& _lambda,
            RooAbsReal& _mu) :
        RooAbsPdf(name,title),
        m_x("x","x",this,_x),
        m_xmin("xmin","xmin",this,_xmin),
        m_lambda("lambda","lambda",this,_lambda),
        m_mu("mu","mu",this,_mu)
{
}

RooInverseGaussian::RooInverseGaussian(const RooInverseGaussian& other, const char* name) :  
    RooAbsPdf(other,name), 
    m_x("x",this,other.m_x),
    m_xmin("xmin",this,other.m_xmin),
    m_lambda("lambda",this,other.m_lambda),
    m_mu("mu",this,other.m_mu)
{ 
} 

// public functions
Double_t RooInverseGaussian::evaluate() const 
{
    Double_t ret = 0.0;
    Double_t xshift = m_x - m_xmin;
    Double_t arg = xshift - m_mu;

    if ( xshift <= 0 ) {
        return ret;
    } else {
        ret = TMath::Sqrt( m_lambda / ( twopi * xshift * xshift * xshift ) ) * 
              TMath::Exp( - m_lambda * arg * arg / ( 2 * m_mu * m_mu * xshift ) );
        return ret;
    }
} 

Int_t RooInverseGaussian::getAnalyticalIntegral(RooArgSet& allVars,
        RooArgSet& analVars, const char* /*rangeName*/) const
{
    if (matchArgs(allVars, analVars, m_x)) return 1;
    return 0;
}

Double_t RooInverseGaussian::analyticalIntegral(Int_t code, const char* rangeName) const
{
    // assert one has to integrate over x
    assert(1 == code);
    Double_t ret = 0.0;

    // integration ranges
    Double_t x0 = m_x.min(rangeName) - m_xmin;
    Double_t x1 = m_x.max(rangeName) - m_xmin;
    Double_t exp2LambdaOverMu = TMath::Exp(2 * m_lambda / m_mu);

    // integrate
    if ( !(x0 > 0.0) ) {
        // just do the cdf
        Double_t x1OverMu = x1 / m_mu ;
        Double_t sqrtLambdaOverX1 = TMath::Sqrt(m_lambda / x1);
        ret = 0.5 * ( 1.0 + RooMath::erf( sqrtLambdaOverX1 * ( x1OverMu - 1.0 ) ) 
                + exp2LambdaOverMu * ( 1.0 + RooMath::erf( - sqrtLambdaOverX1 * ( x1OverMu + 1.0 ) ) ) ); 
        return ret;
    } else {
        // two cdfs
        Double_t x0OverMu = x0 / m_mu ;
        Double_t x1OverMu = x1 / m_mu ;
        Double_t sqrtLambdaOverX0 = TMath::Sqrt(m_lambda / x0);
        Double_t sqrtLambdaOverX1 = TMath::Sqrt(m_lambda / x1);
        ret = 0.5 * ( RooMath::erf( sqrtLambdaOverX1 * ( x1OverMu - 1.0 ) ) - RooMath::erf( sqrtLambdaOverX0 * ( x0OverMu - 1.0 ) )
            + exp2LambdaOverMu * ( RooMath::erf( - sqrtLambdaOverX1 * ( x1OverMu + 1.0 ) ) - RooMath::erf( - sqrtLambdaOverX0 * ( x0OverMu + 1.0 ) ) ) );
        return ret;
    }
}

Int_t RooInverseGaussian::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
    if (matchArgs(directVars,generateVars,m_x)) return 1 ;  
    return 0 ;
}

void RooInverseGaussian::generateEvent(Int_t code)
{
    assert(code==1);
    Double_t xgen, ggen, nu, y, z;
    if(code==1){
        while(1) { 
            nu = RooRandom::gaussian();
            y = nu*nu;
            xgen = m_mu * ( 1.0 + ( m_mu * y - TMath::Sqrt(4.0 * m_mu * m_lambda * y + m_mu * m_mu * y * y) ) / (2.0 * m_lambda) ); 
            z = RooRandom::uniform();
            if ( z > m_mu / ( m_mu + xgen ) ) {
                xgen = m_mu * m_mu / xgen; 
            }
            xgen += m_xmin;
            // ensure range
            if ( xgen < m_x.max() && xgen > m_x.min()) {
                m_x = xgen;
                break;
            }
        }
    } else {
        std::cout << "error in RooInverseGaussian generateEvent"<< std::endl;
    }
    return;
}
