#!/usr/bin/env python3
"""Generate reproducible capture/transport pairs, pilots, thermal runs and cutoff checks."""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path


def prepare(destination: Path, phase: str) -> list[dict]:
    """Write configurations and a manifest; never launch expensive simulations."""
    template=(Path(__file__).resolve().parents[1]/'examples/quickstart.cfg').read_text()
    masses=[.01,.03,.05,.1,.2,.3,.5,1.]
    points={(m,-36,1100) for m in masses}
    points|={(.01,s,1100) for s in [-38,-36,-34,-32,-30]}
    points|={(.1,s,1100) for s in [-38,-34,-30]}
    if phase=='cutoff': points={(.01,-36,r) for r in [1100,3000,10000]}
    if phase=='thermal': points={(m,-36,1100) for m in [2.,3.,4.,10.]}
    destination.mkdir(parents=True,exist_ok=True); manifest=[]
    for i,(mass,exponent,outer) in enumerate(sorted(points)):
        label=f'm{mass:g}_s{exponent}_r{outer}'
        repeats=3 if phase=='production' and (mass,exponent) in {(.01,-36),(.1,-36),(1.,-36),(.01,-38),(.01,-30)} else 1
        for seed_index in range(repeats):
            for kind in (['thermal'] if phase=='thermal' else ['capture','transport']):
                samples=500 if phase!='production' else (100000 if mass in [.01,.1,1.] else 50000)
                # Fixed injection exposure requires pilot-informed tuning at small capture probabilities.
                if kind=='capture': samples=100000 if phase!='production' else 1000000
                seed=20260918+100*i+2*seed_index+(kind=='capture')
                folder=(destination/label/f'seed{seed_index}').resolve()
                replacements={'ID':f'"{label}_{kind}_{seed_index}"','run_mode':'"Capture"' if kind=='capture' else '"Parameter point"',
                              'sample_size':str(samples),'max_trajectories':'0','DM_mass':f'{mass:.12g}',
                              'DM_cross_section_nucleon':f'1.0e{exponent}','output_dir':json.dumps(str(folder)+'/'),
                              'max_trajectory_wall_time_sec':'0.0','interpolation_points':'0'}
                text='// Generated transport configuration; sample_size meaning is specified by workflow.\n'+template[template.index('ID ='):]
                for key,value in replacements.items(): text=re.sub(rf'(?m)^{key}\s*=.*?;',f'{key} = {value};',text)
                text+=f'\nfixed_seed = {seed};\nouter_removal_radius_rsun = {float(outer)};\nproduction_mode = {str(kind!="thermal").lower()};\nthermal_validation_mode = {str(kind=="thermal").lower()};\n'
                if kind=='thermal': text+='maximum_number_of_scatterings = 10000;\n'
                cfg=destination/f'{label}_{kind}_seed{seed_index}.cfg'; cfg.write_text(text)
                manifest.append({'config':str(cfg.resolve()),'kind':kind,'mass_GeV':mass,'sigma_cm2':10.**exponent,'R_remove_rsun':outer,'samples':samples,'seed':seed})
    (destination/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('destination',type=Path)
    parser.add_argument('--phase',choices=['pilot','production','thermal','cutoff'],default='pilot')
    args=parser.parse_args(); print(f'Wrote {len(prepare(args.destination,args.phase))} configs; no jobs launched.')
