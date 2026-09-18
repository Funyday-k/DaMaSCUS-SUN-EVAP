#!/usr/bin/env python3
"""Collisionless, spherical-average unbound SHM reference and finite-volume cross term.

Liouville mapping gives n(r)/n_infinity = integral f_sun(u)*sqrt(1+vesc(r)^2/u^2) du.
This transparent-Sun reference omits capture/attenuation; it is not reconstructed
from the solar-intersecting transit sample, and is not a full sky-background model.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import erf
from scipy.integrate import trapezoid
from analyze_point import R_SUN_CM, require_accepted, physical_config
GM_SUN_CM3_S2=1.32712440018e26


def focused_density(radius_cm: np.ndarray, mass_GeV: float, rho_GeV_cm3: float,
                    v0_kms: float, vobs_kms: float, vgal_kms: float) -> np.ndarray:
    """Exterior monopole unbound density [cm^-3] for a truncated shifted Maxwellian."""
    if (not np.all(np.isfinite(radius_cm)) or not np.all(np.isfinite([mass_GeV,rho_GeV_cm3,v0_kms,vobs_kms,vgal_kms]))
        or np.any(radius_cm<R_SUN_CM) or min(mass_GeV,rho_GeV_cm3,v0_kms,vgal_kms)<=0 or not 0<=vobs_kms<vgal_kms):
        raise ValueError('require exterior radii and 0 <= v_observer < v_galactic_escape')
    u=np.linspace(1e-7,vgal_kms+vobs_kms,12000)
    z=vgal_kms/v0_kms; norm=erf(z)-2*z*np.exp(-z*z)/np.sqrt(np.pi)
    if vobs_kms==0:
        f=4/np.sqrt(np.pi)*u*u/v0_kms**3*np.exp(-(u/v0_kms)**2)/norm
    else:
        tail=np.where(u<vgal_kms-vobs_kms,np.exp(-((u+vobs_kms)/v0_kms)**2),np.exp(-z*z))
        f=u/(np.sqrt(np.pi)*v0_kms*vobs_kms*norm)*(np.exp(-((u-vobs_kms)/v0_kms)**2)-tail)
    f=np.maximum(f,0)
    vesc2=2*GM_SUN_CM3_S2/radius_cm/1e10
    return rho_GeV_cm3/mass_GeV*trapezoid(f[None,:]*np.sqrt(1+vesc2[:,None]/u[None,:]**2),u,axis=1)


def analyze(directory: Path, v0: float, vobs: float, vgal: float) -> dict:
    """Integrate exterior halo, cross, and captured n^2 terms over the recorded domain."""
    meta=require_accepted(directory,'complete_captured_transport')
    derived=json.loads((directory/'derived.json').read_text())
    if derived.get('metadata')!=meta:
        raise ValueError('derived profile belongs to different simulation metadata')
    fraction=float(physical_config(directory)['DM_fraction'])
    profile=np.loadtxt(directory/'tables/radial_profile.tsv')
    low,high=profile[:,0]*R_SUN_CM,profile[:,1]*R_SUN_CM
    sel=low>=R_SUN_CM*(1-1e-12); lo,hi=low[sel],high[sel]
    # Eight-point radial quadrature retains halo variation within each constant-captured-density shell.
    nodes,weights=np.polynomial.legendre.leggauss(8)
    radius=((hi-lo)[:,None]*nodes+(hi+lo)[:,None])/2
    halo=np.stack([focused_density(row,meta['m_chi_GeV'],meta['halo_density_GeV_cm3']*fraction,v0,vobs,vgal) for row in radius])
    dvol=4*np.pi*radius**2*(hi-lo)[:,None]/2*weights
    captured=profile[sel,3,None]
    result={'scope':'transparent-Sun collisionless exterior monopole, finite R_remove volume',
            'DM_fraction':fraction, 'halo_v0_kms':v0,'halo_vobs_kms':vobs,'halo_vgal_kms':vgal,
            'integral_halo2_cm3':float((halo*halo*dvol).sum()),
            'integral_cross_cm3':float((2*halo*captured*dvol).sum()),
            'integral_captured2_cm3':float((captured*captured*dvol).sum())}
    (directory/'halo_reference.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('directory',type=Path)
    p.add_argument('--v0-kms',type=float,required=True); p.add_argument('--vobs-kms',type=float,required=True); p.add_argument('--vgal-kms',type=float,required=True)
    args=p.parse_args(); print(json.dumps(analyze(args.directory,args.v0_kms,args.vobs_kms,args.vgal_kms),indent=2))
