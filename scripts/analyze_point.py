#!/usr/bin/env python3
"""Analyze schema-10 complete histories with an independent fixed-injection capture run.

All densities use cm; trajectory moments are supplied in seconds and km^2/s.
No published detector sensitivity is inferred from a continuum flux threshold.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

R_SUN_CM = 6.957e10
AU_CM = 1.495978707e13
C_KM_S = 299792.458
K_B_EV_K = 8.617333262e-5
BLOCKS = 64


def read_json(path: Path) -> dict:
    """Read metadata without modifying the simulation products."""
    return json.loads(path.read_text())


def normalize_rate_grid_metadata(meta: dict) -> None:
    """Require the explicit effective grid used by schema-10 runs."""
    if any(key not in meta for key in ('rate_radius_points','rate_speed_points','rate_max_speed')):
        raise ValueError('missing scattering-rate grid metadata')
    nr, nv, vmax = (meta[key] for key in
                    ('rate_radius_points', 'rate_speed_points', 'rate_max_speed'))
    if (not isinstance(nr, int) or nr < 0 or not isinstance(nv, int) or nv < 0
        or not np.isfinite(vmax) or vmax < 0
        or ((nr >= 2 and nv >= 2) != (vmax > 0))):
        raise ValueError('invalid scattering-rate grid metadata')


def require_accepted(path: Path, workflow: str) -> dict:
    """Require an accepted, explicitly identified production product."""
    meta = read_json(path / 'metadata.json')
    if meta.get('schema_version') != 10 or meta.get('workflow') != workflow:
        raise ValueError(f'{path}: wrong schema or workflow')
    if not meta.get('seed'):
        raise ValueError(f'{path}: missing reproducible seed')
    if meta.get('production_accepted') is not True:
        raise ValueError(f'{path}: production gate failed; prefixes cannot normalize a source')
    if meta.get('N_numerical_failures')!=0 or meta.get('N_computational_failures')!=0:
        raise ValueError(f'{path}: failure counters contradict accepted production')
    normalize_rate_grid_metadata(meta)
    rank_seeds(meta)
    return meta


def read_capture_result(path: Path, *, log: bool = False) -> dict:
    """Read one stdout record from a Slurm log, or its locally archived JSON."""
    if log:
        prefix = 'CAPTURE_RESULT_JSON='
        with path.open() as stream:
            records = [line[len(prefix):] for line in stream if line.startswith(prefix)]
        if len(records) != 1:
            raise ValueError(f'{path}: expected exactly one {prefix} record, found {len(records)}')
        cap = json.loads(records[0])
    else:
        cap = read_json(path)
    validate_capture_result(cap)
    return cap


def validate_capture_result(cap: dict) -> None:
    """Require accepted capture schema 1, explicit numerics and reproducible RNG streams."""
    if (not isinstance(cap,dict) or cap.get('capture_result_schema') != 1
        or cap.get('workflow') != 'fixed_injection_capture'):
        raise ValueError('wrong capture schema or workflow')
    if cap.get('production_accepted') is not True:
        raise ValueError('capture production gate failed')
    if cap.get('N_numerical_failures') != 0 or cap.get('N_computational_failures') != 0:
        raise ValueError('capture failure counters contradict accepted production')
    normalize_rate_grid_metadata(cap)
    rank_seeds(cap)
    validate_capture_summary(cap)


def rank_seeds(meta: dict) -> set[int]:
    """Reconstruct the uint32 MT19937 seeds, including overflow in C++ rank offsets."""
    seed,ranks=meta.get('seed'),meta.get('mpi_ranks')
    if (not isinstance(seed,int) or not 0<seed<2**32
        or not isinstance(ranks,int) or not 0<ranks<=2**32):
        raise ValueError('invalid master seed or MPI rank count')
    return {(seed+1000003*rank) % 2**32 for rank in range(ranks)}


def validate_capture_summary(cap: dict) -> np.ndarray:
    """Check a complete fixed-incident experiment and its jackknife sufficient statistics."""
    for key in ('N_inj','N_capt','requested_samples'):
        if type(cap.get(key)) is not int or cap[key] < 0:
            raise ValueError('invalid capture counts')
    cb=np.asarray(cap['blocks'],dtype=float)
    if (cb.shape!=(BLOCKS,2) or not np.all(np.isfinite(cb)) or np.any(cb<0)
        or np.any(cb!=np.floor(cb)) or np.any(cb[:,1]>cb[:,0])
        or not np.array_equal(cb.sum(axis=0),[cap['N_inj'],cap['N_capt']])):
        raise ValueError('invalid capture block counts')
    if (cap['fixed_injection'] is not True or cap['N_inj']!=cap['requested_samples']):
        raise ValueError('capture normalization must use a complete fixed incident ensemble')
    if cap['N_capt']<2 or not np.isfinite(cap['C_geom_s_inv']) or cap['C_geom_s_inv']<=0:
        raise ValueError('insufficient capture detections or invalid geometric rate')
    if not np.isclose(cap['f_cap'],cap['N_capt']/cap['N_inj'],rtol=1e-12,atol=0):
        raise ValueError('capture fraction does not match counts')
    if (not np.isfinite(cap['C_capture_s_inv']) or not np.isclose(
            cap['C_capture_s_inv'],cap['C_geom_s_inv']*cap['f_cap'],rtol=1e-12,atol=0)):
        raise ValueError('capture rate does not match geometric rate and fraction')
    return cb


def read_transport_blocks(directory: Path, meta: dict) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Validate native shell geometry and complete block counts before normalization."""
    edges=np.asarray(meta['radial_edges_km'],dtype=float)*1e5
    if (edges.ndim!=1 or len(edges)<3 or not np.all(np.isfinite(edges))
        or edges[0]!=0 or np.any(np.diff(edges)<=0)
        or not np.any(np.isclose(edges,R_SUN_CM,rtol=1e-12,atol=0))
        or not np.isclose(edges[-1]/R_SUN_CM,meta['R_remove_rsun'],rtol=1e-12)
        or meta['R_transit_reference_rsun']>meta['R_remove_rsun']):
        raise ValueError('invalid radial edges or missing solar surface')
    bins=len(edges)-1
    data=np.loadtxt(directory/'radial_blocks.tsv',skiprows=1,ndmin=2)
    if data.shape!=(BLOCKS*bins,12) or not np.all(np.isfinite(data)) or np.any(data[:,4:]<0):
        raise ValueError('incomplete/nonfinite/negative radial block data')
    if not (np.array_equal(data[:,0],np.repeat(np.arange(BLOCKS),bins))
            and np.array_equal(data[:,1],np.tile(np.arange(bins),BLOCKS))):
        raise ValueError('radial block indices must be unique and ordered')
    geometry=np.tile(np.column_stack([edges[:-1],edges[1:]]),(BLOCKS,1))
    if not np.allclose(data[:,2:4]*1e5,geometry,rtol=1e-14,atol=0):
        raise ValueError('metadata and radial edges differ')
    inbound=read_incident_inbound(directory,meta,edges)
    if (not np.allclose(data[:,10].reshape(BLOCKS,bins).sum(axis=0),inbound[:,3],rtol=1e-12,atol=1e-6)
        or not np.allclose(data[:,11].reshape(BLOCKS,bins).sum(axis=0),inbound[:,4],rtol=1e-12,atol=1e-2)):
        raise ValueError('incident inbound block moments do not close')
    samples=np.loadtxt(directory/'block_counts.tsv',skiprows=1,ndmin=2)
    if (samples.shape!=(BLOCKS,2) or not np.all(np.isfinite(samples))
        or not np.array_equal(samples[:,0],np.arange(BLOCKS))
        or np.any(samples[:,1]<0) or np.any(samples[:,1]!=np.floor(samples[:,1]))):
        raise ValueError('invalid transport block counts')
    counts=samples[:,1]
    if (counts.sum()<2 or counts.sum()!=meta['N_completed'] or meta['N_completed']!=meta['N_capt']
        or meta['N_completed']!=meta['requested_samples']):
        raise ValueError('incomplete captured ensemble or captured count mismatch')
    if np.any((counts==0)&(data[:,4].reshape(BLOCKS,bins).sum(axis=1)>0)):
        raise ValueError('empty block contains captured residence')
    return data,edges,counts


def read_incident_inbound(directory: Path, meta: dict, edges_cm: np.ndarray) -> np.ndarray:
    """Validate the complete incoming shell moments against their radial geometry."""
    bins=len(edges_cm)-1
    data=np.loadtxt(directory/'incident_inbound.tsv',skiprows=1,ndmin=2)
    if (data.shape!=(bins,5) or not np.all(np.isfinite(data))
        or np.any(data[:,3:]<0) or not np.array_equal(data[:,0],np.arange(bins))
        or not np.allclose(data[:,1:3]*1e5,
                           np.column_stack([edges_cm[:-1],edges_cm[1:]]),rtol=1e-14,atol=0)
        or (meta['N_inj']>0 and data[:,3].sum()<=0)):
        raise ValueError('invalid incident inbound shell moments')
    return data


def containment(edges_cm: np.ndarray, weight: np.ndarray, fraction: float) -> float:
    """Return a containment radius [R_sun] for piecewise constant volume density."""
    total = weight.sum()
    if total <= 0:
        return float('nan')
    cumulative = np.r_[0., np.cumsum(weight)] / total
    j = min(len(weight)-1, int(np.searchsorted(cumulative, fraction, side='right')-1))
    alpha = (fraction-cumulative[j])/(cumulative[j+1]-cumulative[j])
    return float((edges_cm[j]**3 + alpha*(edges_cm[j+1]**3-edges_cm[j]**3))**(1/3)/R_SUN_CM)


def chord_matrix(edges_cm: np.ndarray, psi: np.ndarray, occulted: bool = True,
                 threshold_cm: float = R_SUN_CM) -> np.ndarray:
    """Exact shell chord lengths [cm] on observer rays l>=0, including r>D_sun."""
    b = AU_CM*np.sin(psi)[:, None]
    center = AU_CM*np.cos(psi)[:, None]
    lo = np.maximum(edges_cm[:-1], threshold_cm)[None, :]
    hi = edges_cm[1:][None, :]
    cap = np.full((len(psi), 1), np.inf)
    if occulted:
        on_disk = (psi < np.pi/2) & (AU_CM*np.sin(psi) < R_SUN_CM)
        cap[on_disk, 0] = AU_CM*np.cos(psi[on_disk])-np.sqrt(R_SUN_CM**2-(AU_CM*np.sin(psi[on_disk]))**2)
    def ball(radius: np.ndarray) -> np.ndarray:
        half = np.sqrt(np.maximum(radius*radius-b*b, 0))
        left, right = np.maximum(0, center-half), np.minimum(cap, center+half)
        return np.where(radius>b, np.maximum(0, right-left), 0)
    return np.where(hi>lo, np.maximum(0, ball(hi)-ball(lo)), 0)


def angular_grid() -> np.ndarray:
    """Resolve the solar limb while retaining the entire sky [radians]."""
    limb = np.arcsin(R_SUN_CM/AU_CM)
    return np.unique(np.r_[0., np.geomspace(limb/10000, limb*.999999, 220),
                           limb, np.geomspace(limb*1.000001, .2, 600),
                           np.linspace(.2, np.pi, 500)])


def cumulative_angle(psi: np.ndarray, intensity: np.ndarray) -> np.ndarray:
    """Integrate intensity over solid angle [ph cm^-2 s^-1]."""
    integrand = 2*np.pi*np.sin(psi)*intensity
    return np.r_[0., np.cumsum(np.diff(psi)*(integrand[:-1]+integrand[1:])/2)]


def observables(dt_s: np.ndarray, v2dt: np.ndarray, n: float, capture_s_inv: float,
                edges: np.ndarray, sigma_v: float, psi: np.ndarray,
                visible: np.ndarray, transparent_weight: np.ndarray,
                distant_weight: np.ndarray) -> tuple[dict, np.ndarray, np.ndarray]:
    """Compute complete-history source observables; sigma_v in cm^3/s."""
    if n <= 0 or capture_s_inv < 0:
        raise ValueError('empty transport ensemble or negative capture rate')
    volume = 4*np.pi/3*np.diff(edges**3)
    tau_bin = dt_s/n
    mu = tau_bin/volume
    quadratic = mu*mu*volume
    inside = edges[1:] <= R_SUN_CM*(1+1e-12)
    result: dict[str, float] = {'C_s_inv': float(capture_s_inv)}
    for name, mask in [('in', inside), ('out', ~inside), ('tot', np.ones(len(mu),dtype=bool))]:
        tau, integral = float(tau_bin[mask].sum()), float(quadratic[mask].sum())
        result.update({f'tau_{name}_s': tau, f'I_{name}_s2_cm3': integral,
                       f'V_eff_{name}_cm3': tau*tau/integral if integral>0 else float('nan'),
                       f'A_{name}_cm3': capture_s_inv**2*integral,
                       f'Gamma_{name}_s_inv': .5*sigma_v*capture_s_inv**2*integral})
    result['f_occ_out'] = result['tau_out_s']/result['tau_tot_s']
    result['f_ann_out'] = result['I_out_s2_cm3']/result['I_tot_s2_cm3']
    result['epsilon_ann'] = sigma_v*capture_s_inv*result['I_tot_s2_cm3']
    for name, weights in [('occ', tau_bin), ('ann', quadratic), ('ann_out',quadratic*(~inside))]:
        for f in [.5,.9]:
            result[f'r{int(100*f)}_{name}_rsun'] = containment(edges,weights,f)
    # 1/2 event factor times two photons per gamma-gamma event.
    emissivity = sigma_v*capture_s_inv**2*mu*mu/(4*np.pi)
    intensity = visible@emissivity
    cumulative = cumulative_angle(psi,intensity)
    flux = float(cumulative[-1])
    transparent = float(transparent_weight@emissivity)
    result.update(Phi_gamma_cm2_s=flux, f_vis=flux/transparent if transparent>0 else float('nan'),
                  eta_beyond_1au=float(distant_weight@emissivity)/flux if flux>0 else float('nan'))
    for f in [.68,.9]:
        result[f'theta{int(f*100)}_deg'] = float(np.interp(f*flux,cumulative,psi)*180/np.pi) if flux>0 else float('nan')
    # All-flavor, neutrino+antineutrino, direct-channel source benchmark. No attenuation/flavor response.
    result['Phi_nu_direct_unattenuated_cm2_s'] = 2*result['Gamma_in_s_inv']/(4*np.pi*AU_CM**2)
    result['v2_in_km2_s2'] = float(v2dt[inside].sum()/dt_s[inside].sum()) if dt_s[inside].sum()>0 else float('nan')
    return result, mu, intensity


def serializable(obj: object) -> object:
    """Represent undefined diagnostics by JSON null, never nonstandard NaN."""
    if isinstance(obj, dict):
        return {k:serializable(v) for k,v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [serializable(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    return obj


def analyze(output: Path, capture: dict, sigma_v: float = 3e-26,
            neutrino_flux_requirement: float | None = None, make_plots: bool = True) -> dict:
    """Write rates, profiles, independent-run jackknife errors, and diagnostic histograms."""
    if not np.isfinite(sigma_v) or sigma_v<=0:
        raise ValueError('sigma_v must be finite and positive')
    m = require_accepted(output,'complete_captured_transport')
    c = capture
    validate_capture_result(c)
    for key in ['m_chi_GeV','sigma_SD_cm2','solar_model','halo_model','halo_density_GeV_cm3',
                'R_inj_rsun','R_match_rsun','R_incident_au',
                'rate_radius_points','rate_speed_points','rate_max_speed',
                'rk_position_tolerance_km','rk_velocity_tolerance_km_s','rk_phase_tolerance',
                'max_optical_depth_step','optical_depth_relative_tolerance']:
        if m[key] != c[key]:
            raise ValueError(f'capture/transport mismatch: {key}')
    if not isinstance(m.get('physical_config'),dict) or not m['physical_config']:
        raise ValueError('transport metadata is missing physical configuration')
    if m['physical_config']!=c.get('physical_config'):
        raise ValueError('capture/transport physical configuration mismatch')
    if rank_seeds(m) & rank_seeds(c):
        raise ValueError('independent capture and transport require disjoint MPI RNG seeds')
    cap = c
    cb = validate_capture_summary(cap)
    data,edges,counts = read_transport_blocks(output,m)
    bins=len(edges)-1
    dt=data[:,4].reshape(BLOCKS,bins); v2=data[:,5].reshape(BLOCKS,bins)
    summary=np.genfromtxt(output/'trajectory_summary.tsv',names=True,dtype=None,encoding='utf8',ndmin=1)
    if len(summary)!=counts.sum() or not np.all(np.isin(summary['termination_reason'],['physical_escape','outer_orbit_removed'])):
        raise ValueError('incomplete captured trajectory summary')
    for key in ['t_end_s','tau_in_s','tau_out_s']:
        if not np.all(np.isfinite(summary[key])) or np.any(summary[key]<0):
            raise ValueError('invalid trajectory residence time')
    block=summary['block_id']
    if (np.any(block!=np.floor(block)) or np.any(block<0) or np.any(block>=BLOCKS)
        or not np.all(summary['seed']==m['seed']) or np.any(summary['rank']<0)
        or np.any(summary['rank']>=m['mpi_ranks'])
        or len(set(zip(summary['rank'],summary['trajectory_id'])))!=len(summary)):
        raise ValueError('invalid or duplicate trajectory identifiers')
    block=block.astype(int)
    if not np.array_equal(np.bincount(block,minlength=BLOCKS),counts):
        raise ValueError('trajectory block counts differ')
    if not np.allclose(summary['tau_in_s']+summary['tau_out_s'],summary['t_end_s'],rtol=1e-8,atol=1e-6):
        raise ValueError('per-trajectory residence accounting failed')
    for field,mask in [('tau_in_s',edges[1:]<=R_SUN_CM*(1+1e-12)),
                       ('tau_out_s',edges[1:]>R_SUN_CM*(1+1e-12))]:
        if not np.allclose(dt[:,mask].sum(axis=1),np.bincount(block,weights=summary[field],minlength=BLOCKS),rtol=1e-8,atol=1e-6):
            raise ValueError('per-block histogram residence accounting failed')
    removed=np.bincount(block,weights=(summary['termination_reason']=='outer_orbit_removed'),minlength=BLOCKS)
    if removed.sum()!=m['N_outer_removed']:
        raise ValueError('removed count differs from trajectory records')
    class_blocks=np.zeros((5,BLOCKS,bins))
    rows=np.loadtxt(output/'orbit_class_blocks.tsv',skiprows=1,ndmin=2)
    if (rows.ndim!=2 or rows.shape[1]!=4 or not np.all(np.isfinite(rows))
        or np.any(rows<0) or np.any(rows[:,:3]!=np.floor(rows[:,:3]))
        or np.any(rows[:,:3]>=np.array([5,BLOCKS,bins]))):
        raise ValueError('invalid aphelion class blocks')
    indices=rows[:,:3].astype(int)
    if len(np.unique(indices,axis=0))!=len(rows):
        raise ValueError('duplicate aphelion class block')
    np.add.at(class_blocks,tuple(indices.T),rows[:,3])
    if not np.allclose(class_blocks.sum(axis=0),dt,rtol=1e-10,atol=1e-8):
        raise ValueError('aphelion class occupation does not close per block')
    classes=class_blocks.sum(axis=1)
    psi=angular_grid(); visible=chord_matrix(edges,psi)
    trapezoid=np.r_[np.diff(psi),0]+np.r_[0,np.diff(psi)]
    solid=trapezoid*np.pi*np.sin(psi)
    trans_weight=solid@chord_matrix(edges,psi,False)
    distant_weight=solid@chord_matrix(edges,psi,True,AU_CM)
    total_dt,total_v2,n=dt.sum(axis=0),v2.sum(axis=0),counts.sum()
    C=cap['C_geom_s_inv']*cap['N_capt']/cap['N_inj']
    def compute(x: np.ndarray,y: np.ndarray,z: float,rate: float) -> tuple:
        return observables(x,y,z,rate,edges,sigma_v,psi,visible,trans_weight,distant_weight)
    point,mu,intensity=compute(total_dt,total_v2,n,C)
    def add_diagnostics(values: dict, fraction: float, removed_fraction: float) -> None:
        values['C_over_C_geom']=fraction
        values['f_removed']=removed_fraction
        values['T2_in_K']=m['m_chi_GeV']*1e9*values['v2_in_km2_s2']/(3*K_B_EV_K*C_KM_S**2)
    add_diagnostics(point,cap['N_capt']/cap['N_inj'],removed.sum()/n)
    keys=list(point)
    transport_reps=[]; capture_reps=[]
    for k in range(BLOCKS):
        if n-counts[k]<=0 or cap['N_inj']-cb[k,0]<=0:
            raise ValueError('too few histories for delete-block uncertainty')
        values=compute(total_dt-dt[k],total_v2-v2[k],n-counts[k],C)[0]
        add_diagnostics(values,cap['N_capt']/cap['N_inj'],(removed.sum()-removed[k])/(n-counts[k]))
        transport_reps.append([values[j] for j in keys])
        rate=cap['C_geom_s_inv']*(cap['N_capt']-cb[k,1])/(cap['N_inj']-cb[k,0])
        values=compute(total_dt,total_v2,n,rate)[0]
        add_diagnostics(values,rate/cap['C_geom_s_inv'],removed.sum()/n)
        capture_reps.append([values[j] for j in keys])
    def variance(x: list) -> np.ndarray:
        a=np.asarray(x); return (BLOCKS-1)/BLOCKS*((a-a.mean(axis=0))**2).sum(axis=0)
    errors=np.sqrt(variance(transport_reps)+variance(capture_reps))
    transport_array=np.asarray(transport_reps,dtype=float)
    capture_array=np.asarray(capture_reps,dtype=float)
    central_array=np.asarray([point[key] for key in keys],dtype=float)
    # Delete-block bias is evaluated independently for transport and capture,
    # then added. In particular, I=integral(mu^2 dV) has positive plug-in bias.
    bias=(BLOCKS-1)*(transport_array.mean(axis=0)-central_array)
    bias+=(BLOCKS-1)*(capture_array.mean(axis=0)-central_array)
    corrected=central_array-bias
    result={'analysis_version':2,'physical_config':m['physical_config'],
            'metadata':m,'capture_metadata':c,'sigma_v_cm3_s':sigma_v,'central':point,
            'jackknife_se':dict(zip(keys,errors)),
            'jackknife_bias':dict(zip(keys,bias)),
            'jackknife_bias_corrected':dict(zip(keys,corrected)),
            'uncertainty':'independent capture and transport delete-block variances added; each replicate recomputes all nonlinear observables',
            'neutrino_model':'direct nu+antinu all-flavor source benchmark, unattenuated, no detector likelihood'}
    if neutrino_flux_requirement is not None:
        result['Gamma_nu_required_s_inv']=4*np.pi*AU_CM**2*neutrino_flux_requirement/2
        result['neutrino_flux_requirement_cm2_s']=neutrino_flux_requirement
    tables=output/'tables'; tables.mkdir(exist_ok=True)
    volume=4*np.pi/3*np.diff(edges**3)
    np.savetxt(tables/'radial_profile.tsv',np.column_stack([edges[:-1]/R_SUN_CM,edges[1:]/R_SUN_CM,mu,C*mu,total_dt/n,mu*mu*volume]),delimiter='\t',header='r_low_rsun r_high_rsun mu_s_cm3 n_cm3 tau_bin_s I_bin_s2_cm3')
    np.savetxt(tables/'gamma_profile.tsv',np.column_stack([np.degrees(psi),intensity,cumulative_angle(psi,intensity)]),delimiter='\t',header='psi_deg intensity_cm2_s_sr cumulative_flux_cm2_s')
    # Additive attribution includes cross-class annihilation terms: integral mu_c*mu_total.
    class_mu=classes/n/volume
    out=edges[:-1]>=R_SUN_CM*(1-1e-12)
    attribution=[{'class':i,'tau_out_s':float((classes[i,out]/n).sum()),
                  'I_out_attributed_s2_cm3':float((class_mu[i,out]*mu[out]*volume[out]).sum())} for i in range(5)]
    result['aphelion_classes']=attribution
    for name,column in [('t_end',summary['t_end_s']),('r_ap',summary['max_aphelion_km']/R_SUN_CM*1e5),('n_scatter',summary['n_scatter'])]:
        x=np.asarray(column,dtype=float); x=x[np.isfinite(x)&(x>0)]
        if len(x):
            h,e=np.histogram(x,bins=np.geomspace(x.min()*.999,x.max()*1.001,41))
            np.savetxt(output/f'hist_{name}.tsv',np.column_stack([e[:-1],e[1:],h]),delimiter='\t',header='lower upper count')
    source={'mass_GeV':m['m_chi_GeV'],'sigma_v_cm3_s':sigma_v,
            'Gamma_in_s_inv':point['Gamma_in_s_inv'],'Gamma_out_s_inv':point['Gamma_out_s_inv'],
            'Gamma_total_s_inv':point['Gamma_tot_s_inv'],'channel':'unit-branching reference; select gamma-gamma or nu-antinu separately'}
    (output/'annihilation_source.json').write_text(json.dumps(source,indent=2)+'\n')
    (output/'derived.json').write_text(json.dumps(serializable(result),indent=2,allow_nan=False)+'\n')
    if make_plots:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,2,figsize=(9,7))
        centers=np.sqrt(edges[:-1]*edges[1:])/R_SUN_CM
        dln=np.log(edges[1:]/np.maximum(edges[:-1],edges[1:]*1e-12))
        axes[0,0].loglog(centers,total_dt/n/point['tau_tot_s']/dln)
        axes[0,0].set(xlabel=r'$r/R_\odot$',ylabel=r'$dP_{occ}/d\ln r$')
        axes[0,1].loglog(centers,mu*mu*volume/point['I_tot_s2_cm3']/dln)
        axes[0,1].set(xlabel=r'$r/R_\odot$',ylabel=r'$I^{-1}dI/d\ln r$')
        axes[1,0].loglog(np.degrees(psi[1:]),intensity[1:])
        axes[1,0].set(xlabel='Angle [deg]',ylabel=r'Intensity [cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        axes[1,1].semilogx(np.degrees(psi[1:]),cumulative_angle(psi,intensity)[1:]/point['Phi_gamma_cm2_s'])
        axes[1,1].set(xlabel='Angle [deg]',ylabel='Visible cumulative fraction',ylim=(0,1.02))
        fig.tight_layout(); dest=output/'figures'; dest.mkdir(exist_ok=True); fig.savefig(dest/'point_summary.pdf'); plt.close(fig)
    return result


def main() -> None:
    """CLI; rates are conditional on the supplied annihilation coefficient."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    capture_group=parser.add_mutually_exclusive_group(required=True)
    capture_group.add_argument('--capture-log',type=Path)
    capture_group.add_argument('--capture-json',type=Path)
    parser.add_argument('--sigma-v',type=float,default=3e-26)
    parser.add_argument('--neutrino-flux-requirement',type=float)
    parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args()
    if not np.isfinite(args.sigma_v) or args.sigma_v<=0:
        parser.error('--sigma-v must be finite and positive')
    if args.neutrino_flux_requirement is not None and (not np.isfinite(args.neutrino_flux_requirement) or args.neutrino_flux_requirement<=0):
        parser.error('--neutrino-flux-requirement must be finite and positive')
    capture=read_capture_result(args.capture_log or args.capture_json,log=args.capture_log is not None)
    result=analyze(args.output,capture,args.sigma_v,args.neutrino_flux_requirement,not args.no_plots)
    print(json.dumps(serializable(result['central']),indent=2))

if __name__=='__main__':
    main()
