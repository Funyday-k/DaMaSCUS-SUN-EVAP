#!/usr/bin/env python3
"""Exercise stdout-only capture, reduced transport, and local diagnostics with real MPI runs."""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def main() -> None:
    """Run small deterministic configurations; time budgets are in wall-clock seconds."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--program',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True)
    parser.add_argument('--mpiexec')
    parser.add_argument('--numproc-flag',default='-np')
    parser.add_argument('--ranks',type=int,default=1)
    parser.add_argument('--analyze',action='store_true')
    args=parser.parse_args()
    repo=Path(__file__).resolve().parents[1]
    template=(repo/'examples/quickstart.cfg').read_text()
    launcher=[args.mpiexec,args.numproc_flag,str(args.ranks)] if args.mpiexec else []
    science={'metadata.json','radial_blocks.tsv','block_counts.tsv','trajectory_summary.tsv',
             'orbit_class_blocks.tsv','incident_inbound.tsv','termination_counts.tsv','solar_reference.tsv','snapshot'}
    with tempfile.TemporaryDirectory(prefix='damascus-production-contract-') as tmp:
        root=Path(tmp)

        def run(name: str, **changes: str) -> tuple[subprocess.CompletedProcess,Path,Path]:
            """Write an external cfg and launch it, preserving stdout for the analyzer."""
            output=root/name
            settings={'run_mode':'"Capture"','sample_size':'128','max_trajectories':'0',
                      'DM_cross_section_nucleon':'1.0e-32','fixed_seed':'20260921',
                      'production_mode':'true','max_trajectory_wall_time_sec':'0.0',
                      'output_dir':json.dumps(str(output)+'/'),'snapshot_enabled':'true','snapshot_interval':'1',
                      # Production must suppress these even if an old cfg enables them.
                      'trajectory_summary_enabled':'true','trajectory_events_enabled':'true',
                      'trajectory_trace_rate':'1.0',**changes}
            text=template
            for key,value in settings.items():
                pattern=rf'(?m)^{key}\s*=.*?;'
                if re.search(pattern,text): text=re.sub(pattern,f'{key} = {value};',text)
                else: text+=f'\n{key} = {value};\n'
            config=root/f'{name}.cfg'; config.write_text(text)
            result=subprocess.run(launcher+[str(args.program),str(config)],cwd=root,text=True,
                                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
            (root/f'{name}.out').write_text(result.stdout)
            return result,output,config

        def capture_record(result: subprocess.CompletedProcess) -> dict:
            """Check global integer block sums and the normalization printed by rank zero."""
            lines=[line.split('=',1)[1] for line in result.stdout.splitlines() if line.startswith('CAPTURE_RESULT_JSON=')]
            assert len(lines)==1,result.stdout
            cap=json.loads(lines[0])
            assert cap['capture_result_schema']==1 and cap['mpi_ranks']==args.ranks
            assert len(cap['blocks'])==64
            assert cap['N_inj']==sum(row[0] for row in cap['blocks'])
            assert cap['N_capt']==sum(row[1] for row in cap['blocks'])
            assert math.isclose(cap['f_cap'],cap['N_capt']/cap['N_inj'])
            assert math.isclose(cap['C_capture_s_inv'],cap['C_geom_s_inv']*cap['f_cap'])
            return cap

        result,output,_=run('capture')
        assert result.returncode==0,result.stdout
        cap=capture_record(result)
        assert cap['production_accepted'] is True and cap['N_inj']==128
        assert cap['N_numerical_failures']==cap['N_computational_failures']==0
        assert not output.exists(), 'Capture created its output directory'

        result,output,_=run('failed_capture',max_trajectory_wall_time_sec='1e-12')
        assert result.returncode==2,result.stdout
        failed=capture_record(result)
        assert failed['production_accepted'] is False and failed['N_computational_failures']>0
        assert failed['N_inj']<=args.ranks
        assert not output.exists(), 'Failed capture created output files'

        # Seed one stale local run's artifacts to verify production output isolation on rerun.
        product=root/'transport/results_-2.000000_-32.000000'
        product.mkdir(parents=True)
        for name in ['input.cfg','run_metadata.json','diagnostic_trajectory_summary.tsv','trajectory_events.tsv',
                     'invalid_trajectories.tsv','bincount.txt','capture_summary.json','evaporation_times.txt',
                     'residence_jackknife_blocks.tsv']:
            (product/name).write_text('stale local diagnostic')
        result,_,_=run('transport',run_mode='"Parameter point"',sample_size='64',fixed_seed='20260922')
        assert result.returncode==0,result.stdout
        assert {path.name for path in product.iterdir()}==science
        meta=json.loads((product/'metadata.json').read_text())
        assert meta['schema_version']==10 and meta['production_accepted'] is True
        assert meta['N_completed']==meta['N_capt']==64
        assert meta['maximum_number_of_scatterings']>0 and meta['max_trajectory_wall_time_sec']==0
        assert meta['production_mode'] is True and meta['thermal_validation_mode'] is False
        assert meta['restart_supported'] is False
        assert (product/'snapshot').is_dir()

        result,output,_=run('failed_transport',run_mode='"Parameter point"',sample_size='4',
                            max_trajectory_wall_time_sec='1e-12')
        assert result.returncode==2,result.stdout
        failed_product=output/'results_-2.000000_-32.000000'
        assert {path.name for path in failed_product.iterdir()}==science
        assert json.loads((failed_product/'metadata.json').read_text())['production_accepted'] is False

        result,output,_=run('local',run_mode='"Parameter point"',sample_size='4',
                            fixed_seed='20260922',production_mode='false',snapshot_enabled='false')
        assert result.returncode==0,result.stdout
        local=output/'results_-2.000000_-32.000000'
        for name in ['run_metadata.json','diagnostic_trajectory_summary.tsv','trajectory_events.tsv','invalid_trajectories.tsv']:
            assert (local/name).is_file(),name
        rows=list(csv.DictReader(io.StringIO((local/'diagnostic_trajectory_summary.tsv').read_text()),delimiter='\t'))
        assert rows and any(row['rng_state_before_simulation'] for row in rows)
        assert len((local/'trajectory_events.tsv').read_text().splitlines())>1

        result,output,config=run('local_invalid',run_mode='"Parameter point"',sample_size='1',
                                 max_trajectories='1',maximum_number_of_scatterings='1',
                                 production_mode='false',snapshot_enabled='false')
        assert result.returncode==0,result.stdout
        ledger=output/'results_-2.000000_-32.000000/invalid_trajectories.tsv'
        rows=list(csv.DictReader((line for line in ledger.read_text().splitlines() if not line.startswith('#')),delimiter='\t'))
        assert rows and rows[0]['rng_state_before_simulation']
        replay=subprocess.run([str(args.replay),str(config),str(ledger),rows[0]['rank'],rows[0]['trajectory_id']],
                              cwd=root,capture_output=True,text=True,timeout=20)
        assert replay.returncode==0,replay.stdout+replay.stderr
        assert 'original_termination_reason=max_scatterings' in replay.stdout
        assert 'replay_scatterings=1' in replay.stdout

        if args.analyze:
            result=subprocess.run([sys.executable,str(repo/'scripts/analyze_point.py'),str(product),
                                   '--capture-log',str(root/'capture.out'),'--no-plots'],
                                  cwd=root,capture_output=True,text=True,timeout=30)
            assert result.returncode==0,result.stderr
            derived=json.loads((product/'derived.json').read_text())
            assert math.isclose(derived['central']['C_s_inv'],cap['C_capture_s_inv'])
            for key in ['C_s_inv','Gamma_tot_s_inv','I_out_s2_cm3']:
                assert math.isfinite(derived['central'][key])
                assert math.isfinite(derived['jackknife_se'][key])
        print(f'{args.ranks} rank(s): capture, failure, transport, local replay and analysis contracts passed')


if __name__=='__main__':
    main()
