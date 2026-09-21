#!/usr/bin/env python3
"""Collect compatible accepted derived.json files with independent seed uncertainties."""
from __future__ import annotations
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from analyze_point import normalize_rate_grid_metadata, rank_seeds


def comparison_signature(record: dict) -> str:
    """Keep every source/model setting except the three explicitly scanned axes."""
    if record.get('analysis_version')!=2 or not record.get('physical_config') or not record.get('solar_reference_sha256'):
        raise ValueError('missing analysis provenance; rerun analyze_point.py')
    meta=record['metadata']
    normalize_rate_grid_metadata(meta)
    physical={k:v for k,v in record['physical_config'].items()
              if k not in {'DM_mass','DM_cross_section_nucleon'}}
    settings={key:meta[key] for key in ['source_sha256','solar_model','halo_model','halo_density_GeV_cm3',
              'R_inj_rsun','R_match_rsun','R_incident_au','interpolation_points',
              'rate_radius_points','rate_speed_points','rate_max_speed',
              'rk_position_tolerance_km','rk_velocity_tolerance_km_s','rk_phase_tolerance',
              'max_optical_depth_step','optical_depth_relative_tolerance']}
    settings.update(physical=physical,solar_reference=record['solar_reference_sha256'],
                    sigma_v=record['sigma_v_cm3_s'],analysis=record['analysis_version'])
    return json.dumps(settings,sort_keys=True,allow_nan=False)


def summarize(files: list[Path], destination: Path) -> None:
    """Compare mass/cross-section/cutoff axes only; require independent pairs at each point."""
    if not files:
        raise ValueError('no derived products supplied')
    groups=defaultdict(list)
    signature=None
    for file in files:
        record=json.loads(file.read_text()); m=record['metadata']; c=record['capture_metadata']
        for meta,workflow in [(m,'complete_captured_transport'),(c,'fixed_injection_capture')]:
            if (meta.get('production_accepted') is not True or meta.get('schema_version')!=8
                or meta.get('workflow')!=workflow or meta.get('N_numerical_failures')!=0
                or meta.get('N_computational_failures')!=0):
                raise ValueError(f'{file}: rejected production')
        current=comparison_signature(record)
        if signature is not None and current!=signature:
            raise ValueError('incompatible analysis coefficient, physical model, or source provenance')
        signature=current
        if rank_seeds(m) & rank_seeds(c):
            raise ValueError('capture/transport MPI RNG streams overlap')
        groups[(m['m_chi_GeV'],m['sigma_SD_cm2'],m['R_remove_rsun'])].append(record)
    rows=[]
    for (mass,sigma,outer),group in sorted(groups.items()):
        used=set()
        for record in group:
            streams=rank_seeds(record['metadata']) | rank_seeds(record['capture_metadata'])
            if used & streams:
                raise ValueError('seed repeats or shared capture normalization across alleged independent runs')
            used.update(streams)
        row={'mass_GeV':mass,'sigma_cm2':sigma,'R_remove_rsun':outer,'seeds':len(group),
             'sigma_v_cm3_s':group[0]['sigma_v_cm3_s']}
        for key in group[0]['central']:
            vals=[g['central'][key] for g in group]
            if any(v is None for v in vals):
                continue
            errors=[g['jackknife_se'].get(key) for g in group]
            if (not np.all(np.isfinite(vals)) or any(v is None for v in errors)
                or not np.all(np.isfinite(errors)) or min(errors)<0):
                raise ValueError(f'missing/nonfinite uncertainty for {key}; rerun point analysis')
            row[key]=float(np.mean(vals))
            seed=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
            row[key+'_se']=max(max(errors),seed)
            biases=[g['jackknife_bias'].get(key) for g in group]
            if any(v is None or not np.isfinite(v) for v in biases):
                raise ValueError(f'missing/nonfinite jackknife bias for {key}')
            row[key+'_jackknife_bias']=float(np.mean(biases))
        rows.append(row)
    slopes=[]
    for mass,outer in sorted({(r['mass_GeV'],r['R_remove_rsun']) for r in rows}):
        scan=sorted((r for r in rows if r['mass_GeV']==mass and r['R_remove_rsun']==outer),key=lambda r:r['sigma_cm2'])
        for a,b in zip(scan,scan[1:]):
            rec={'mass_GeV':mass,'R_remove_rsun':outer,'sigma_low':a['sigma_cm2'],'sigma_high':b['sigma_cm2']}
            for key in ['C_s_inv','tau_out_s','I_out_s2_cm3','V_eff_out_cm3','A_out_cm3']:
                x,y=a.get(key),b.get(key)
                rec['s_'+key]=(float(np.log(y/x)/np.log(b['sigma_cm2']/a['sigma_cm2']))
                               if x is not None and y is not None and x>0 and y>0 else None)
            terms=[rec['s_A_out_cm3'],rec['s_C_s_inv'],rec['s_I_out_s2_cm3']]
            rec['s_A_minus_2s_C_minus_s_I']=terms[0]-2*terms[1]-terms[2] if all(t is not None for t in terms) else None
            slopes.append(rec)
    destination.mkdir(parents=True,exist_ok=True)
    for name,records in [('points.csv',rows),('slopes.csv',slopes)]:
        # Write an empty file as well, so an old slope table cannot survive a one-point rerun.
        with (destination/name).open('w',newline='') as stream:
            if records:
                fields=sorted(set().union(*(r.keys() for r in records)))
                writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); writer.writerows(records)
    (destination/'manifest.json').write_text(json.dumps({'inputs':[str(f.resolve()) for f in files],
        'comparison_settings':json.loads(signature),
        'seed_error':'max(largest within-run JK SE, independent-pair sample SD)'},indent=2)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('files',nargs='+',type=Path); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); summarize(args.files,args.output)
