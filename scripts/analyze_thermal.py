#!/usr/bin/env python3
"""Interior shape validation against the isothermal SD-proton energy-balance reference.

Reads the potential and target profiles exported by the same C++ solar model.
This workflow never produces absolute occupation, annihilation rates, or fluxes.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import brentq
from analyze_point import R_SUN_CM, K_B_EV_K, C_KM_S, BLOCKS, read_json, read_transport_blocks, serializable


def analyze(directory: Path) -> dict:
    """Return D_TV and T2/Tchi for one explicitly labelled thermal sample."""
    meta=read_json(directory/'metadata.json')
    if (meta.get('schema_version')!=8 or meta.get('workflow')!='thermal_shape_validation'
        or meta.get('N_numerical_failures')!=0 or meta.get('production_accepted') is not False):
        raise ValueError('require a numerical-failure-free thermal shape workflow')
    table=np.loadtxt(directory/'solar_reference.tsv',skiprows=1)
    if (table.ndim!=2 or table.shape[1]!=4 or not np.all(np.isfinite(table))
        or np.any(np.diff(table[:,0])<=0) or np.any(table[:,1]<=0) or np.any(table[:,2:]<0)):
        raise ValueError('invalid solar reference table')
    radius,T,nH,phi=table.T
    mass=meta['m_chi_GeV']*1e9  # rest energy [eV]
    proton=0.93827208816e9
    def balance(Tchi: float) -> float:
        boltz=np.exp(-mass*phi/(C_KM_S**2*K_B_EV_K*Tchi))
        speed=np.sqrt(K_B_EV_K*(Tchi/mass+T/proton))
        return float(trapezoid(nH*speed*(T-Tchi)*boltz*radius**2,radius))
    Tchi=brentq(balance,.3*T.min(),1.5*T.max(),xtol=1e-6,rtol=1e-10)
    radial=radius**2*np.exp(-mass*phi/(C_KM_S**2*K_B_EV_K*Tchi))
    cdf=cumulative_trapezoid(radial,radius,initial=0); cdf/=cdf[-1]
    data,edges,counts=read_transport_blocks(directory,meta)
    inside=edges[1:]<=R_SUN_CM*(1+1e-12)
    hist=data[:,4].reshape(BLOCKS,-1)[:,inside]; v2=data[:,5].reshape(BLOCKS,-1)[:,inside]
    if hist.sum()<=0 or np.any(hist.sum()-hist.sum(axis=1)<=0):
        raise ValueError('insufficient independent interior residence for thermal uncertainty')
    prob=hist.sum(axis=0)/hist.sum()
    expected=np.diff(np.interp(edges[:len(prob)+1],radius,cdf)); expected/=expected.sum()
    def measures(h: np.ndarray,v: np.ndarray) -> np.ndarray:
        return np.array([.5*np.abs(h/h.sum()-expected).sum(),mass*v.sum()/h.sum()/(3*C_KM_S**2*K_B_EV_K*Tchi)])
    central=measures(hist.sum(axis=0),v2.sum(axis=0))
    reps=np.array([measures(hist.sum(axis=0)-h,v2.sum(axis=0)-v) for h,v in zip(hist,v2)])
    se=np.sqrt((BLOCKS-1)/BLOCKS*((reps-reps.mean(axis=0))**2).sum(axis=0))
    result={'m_chi_GeV':meta['m_chi_GeV'],'D_TV':central[0],'T2_over_Tchi':central[1],
            'Tchi_K':Tchi,'jackknife_se_D_TV':se[0],'jackknife_se_T2_over_Tchi':se[1],
            'meaning':'interior shape only; vary scattering budget and independent histories to test relaxation',
            'production_accepted':False}
    (directory/'thermal_validation.json').write_text(json.dumps(serializable(result),indent=2,allow_nan=False)+'\n')
    np.savetxt(directory/'thermal_profiles.tsv',np.column_stack([edges[1:len(prob)+1]/R_SUN_CM,prob,expected,np.cumsum(prob-expected)]),delimiter='\t',header='r_upper_rsun occupation_probability isothermal_probability cumulative_residual')
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('directory',type=Path)
    print(json.dumps(analyze(parser.parse_args().directory),indent=2))
