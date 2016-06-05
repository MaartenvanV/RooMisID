/*****************************************************************************
 * Project: Momentum fraction pdf for mis ID invariant mass correction       *
 *                                                                           *
 * @author: Maarten van Veghel                                               * 
 * @date: 01-2016                                                            *
 *****************************************************************************/

#include "RooFit.h"
#include "Riostream.h" 

#include <cassert>
#include <cmath>
#include "TMath.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooMsgService.h"

#include "RooPowerLawApproxVar.h"

// constructors
RooPowerLawApproxVar::RooPowerLawApproxVar(const char *name, const char *title,
            RooAbsReal& _power,
            const Int_t _index,
            const Double_t _beta_approx,
            const Int_t _n_expos,
            Flags _flags) :
    RooAbsReal(name,title),
    m_index(_index),
    m_power("power","power",this,_power),
    m_beta_approx(_beta_approx),
    m_n_expos(_n_expos),
    m_flags(_flags)
{
}

RooPowerLawApproxVar::RooPowerLawApproxVar(const RooPowerLawApproxVar& other, const char* name) :
    RooAbsReal(other,name),
    m_index(other.m_index),
    m_power("power",this,other.m_power),
    m_beta_approx(other.m_beta_approx),
    m_n_expos(other.m_n_expos),
    m_flags(other.m_flags)
{
}

Double_t RooPowerLawApproxVar::evaluate() const 
{
    // evaluation function
    if ( m_flags & Coeff ) {
        // coefficient of exponent
        Double_t sum = 0.0;
        Double_t beta = m_beta_approx;
        Double_t alph = m_power;
        Int_t idx = m_index;
        for ( Int_t i = 0; i < (Int_t)m_n_expos; i++ ) {
            sum += exp(-alph/TMath::Power(beta,i)) / TMath::Power(beta,i*alph);
        }
        return 1.0 / TMath::Power(beta,idx*alph) / sum ;
    } 
    else if ( m_flags & Slope ) {
        // slope of exponent
        return - ((Double_t)m_power) / TMath::Power(((Double_t)m_beta_approx),((Int_t)m_index)) ;
    } 
    else {
        // must not go here
        assert(1==0);
    }
}

