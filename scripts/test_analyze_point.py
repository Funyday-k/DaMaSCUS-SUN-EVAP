#!/usr/bin/env python3
"""Independent geometry and production-rejection tests for transport analysis."""
import json
import copy
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
from analyze_point import analyze, read_capture_result, population_density, population_annihilation, read_transport_blocks, BLOCKS, K_B_EV_K, C_KM_S
from summarize_transport_scan import summarize
from analyze_point import AU_CM, R_SUN_CM, chord_matrix, containment, require_accepted, angular_grid, cumulative_angle

class GeometryTests(unittest.TestCase):
    def test_equal_populations_include_annihilation_cross_term(self) -> None:
        """Two equal densities give four times either self-pair rate [s^-1]."""
        edges=np.array([0.,R_SUN_CM,2*R_SUN_CM])
        rates=population_annihilation(np.full((2,3),[2.,2.,4.]),edges,3e-26)
        np.testing.assert_allclose(rates[:,3],4*rates[:,0])
        np.testing.assert_allclose(rates[:,2],2*rates[:,0])
        np.testing.assert_allclose(rates[:,3],.5*3e-26*16*(4*np.pi/3*np.diff(edges**3)))

    def test_annihilation_preserves_native_density_variation(self) -> None:
        """Equal-volume shells with n=0,2 have twice the rate of their mean n=1."""
        edges=np.array([0.,1.,np.cbrt(2.)])*R_SUN_CM
        rates=population_annihilation(np.array([[0.,0.,0.],[2.,0.,2.]]),edges,3e-26)
        averaged=population_annihilation(np.array([[1.,0.,1.]]),edges[[0,2]],3e-26)
        np.testing.assert_allclose(rates[:,3].sum(),2*averaged[0,3])

    def test_observer_inside_a_shell_sees_sources_behind_them(self) -> None:
        edges=np.array([0.0,2*AU_CM,3*AU_CM]); psi=np.array([0.0,np.pi/2,np.pi])
        path=chord_matrix(edges,psi,False,2*AU_CM)
        self.assertAlmostEqual(path[2,1]/AU_CM,1.0,places=12)
        self.assertAlmostEqual(path[0,1]/AU_CM,1.0,places=12)
        self.assertAlmostEqual(path[1,1]/AU_CM,np.sqrt(8)-np.sqrt(3),places=12)

    def test_occultation_preserves_near_side_and_removes_far_side(self) -> None:
        edges=np.array([R_SUN_CM,2*R_SUN_CM])
        transparent=chord_matrix(edges,np.array([0.0]),False)
        visible=chord_matrix(edges,np.array([0.0]),True)
        self.assertAlmostEqual(transparent[0,0]/R_SUN_CM,2.0,places=10)
        self.assertAlmostEqual(visible[0,0]/R_SUN_CM,1.0,places=10)

    def test_full_sky_flux_against_volume_geometry_beyond_observer(self) -> None:
        # A constant emissivity shell 2D<r<3D: analytic angular integration
        # over source-centered shells gives (2*pi/D)*integral r*log((r+D)/(r-D)) dr.
        from scipy.integrate import quad
        edges=np.array([2*AU_CM,3*AU_CM]); psi=angular_grid()
        expected=2*np.pi*AU_CM*quad(lambda x:x*np.log((x+1)/(x-1)),2,3,epsabs=1e-11)[0]
        computed=cumulative_angle(psi,chord_matrix(edges,psi,False,2*AU_CM)[:,0])[-1]
        self.assertLess(abs(computed/expected-1),1e-5)

    def test_containment_uses_volume_not_radius(self) -> None:
        result=containment(np.array([0.,R_SUN_CM]),np.array([1.]),.5)
        self.assertAlmostEqual(result,2**(-1/3),places=12)

    def test_unfinished_histories_and_wrong_workflows_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/'metadata.json').write_text(json.dumps({'schema_version':10,'workflow':'complete_captured_transport','production_accepted':False}))
            with self.assertRaises(ValueError): require_accepted(p,'complete_captured_transport')
            (p/'metadata.json').write_text(json.dumps({'schema_version':10,'workflow':'thermal_shape_validation','production_accepted':True}))
            with self.assertRaises(ValueError): require_accepted(p,'complete_captured_transport')

class AnalysisContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name); self.output=self.root/'transport'; self.capture=self.root/'capture.json'
        self.output.mkdir()
        self.meta={'schema_version':10,'production_accepted':True,'workflow':'complete_captured_transport',
            'source_sha256':'fixture-source','seed':101,'mpi_ranks':1,'solar_model':'fixture',
            'halo_model':'SHM','halo_density_GeV_cm3':.4,'m_chi_GeV':.01,'sigma_SD_cm2':1e-32,
            'physical_config':{'DM_mass':.01,'DM_cross_section_nucleon':1e-32,
                               'DM_fraction':1.0,'DM_distribution':'SHM'},
            'R_inj_rsun':2.0,'R_match_rsun':1.0,'R_incident_au':1100,
            'R_transit_reference_rsun':3,'R_remove_rsun':3,
            'interpolation_points':0,'rate_radius_points':0,'rate_speed_points':0,
            'rate_max_speed':0.0,'rate_query_count':100,'rate_fallback_count':0,
            'rate_fallback_fraction':0.0,'rate_max_speed_seen':0.01,
            'rk_position_tolerance_km':1.0,
            'rk_velocity_tolerance_km_s':.001,'rk_phase_tolerance':1e-7,
            'max_optical_depth_step':.05,'optical_depth_relative_tolerance':.01,
            'requested_samples':64,
            'N_completed':64,'N_capt':64,'N_inj':128,'N_outer_removed':32,
            'N_numerical_failures':0,'N_computational_failures':0,
            'radial_edges_km':[0,R_SUN_CM/1e5,2*R_SUN_CM/1e5,3*R_SUN_CM/1e5]}
        self.cmeta={**self.meta,'workflow':'fixed_injection_capture','seed':202,'requested_samples':1280,
                    'N_inj':1280,'N_capt':672,'N_completed':0}
        self.cap={'fixed_injection':True,'N_inj':1280,'N_capt':672,'C_geom_s_inv':1e20,
                  'f_cap':672/1280,'blocks':[[20,10+k%2] for k in range(BLOCKS)]}
        self.cap={**self.cmeta,**self.cap,'capture_result_schema':1,'C_capture_s_inv':1e20*672/1280}
        del self.cap['schema_version']
        self.cmeta=self.cap
        self.write_json(self.output/'metadata.json',self.meta)
        (self.output/'solar_reference.tsv').write_text('transport solar reference')
        self.write_json(self.capture,self.cap)
        blocks=[]; classes=[]; records=[]
        self.tin=1+np.arange(BLOCKS)%3; self.tout=2+np.arange(BLOCKS)%5
        self.speed2=100000+1000*np.arange(BLOCKS)
        for k in range(BLOCKS):
            for b,dt in enumerate([self.tin[k],self.tout[k]/2,self.tout[k]/2]):
                blocks.append([k,b,self.meta['radial_edges_km'][b],self.meta['radial_edges_km'][b+1],dt,dt*self.speed2[k],0,0,0,0,1 if b else 0,100 if b else 0])
                classes.append([0,k,b,dt])
            reason='outer_orbit_removed' if k%2 else 'physical_escape'
            records.append(f'{k+1} 0 101 {k} {reason} {self.tin[k]+self.tout[k]} {self.tin[k]} {self.tout[k]} 2000000 {k+1}')
        np.savetxt(self.output/'radial_blocks.tsv',blocks,header='block bin lo hi dt v2dt transit_dt transit_v2dt post_dt post_v2dt inbound_dt inbound_v2dt')
        for directory,scale in [(self.output,BLOCKS)]:
            inbound=[[b,self.meta['radial_edges_km'][b],self.meta['radial_edges_km'][b+1],
                      scale if b else 0,100*scale if b else 0] for b in range(3)]
            np.savetxt(directory/'incident_inbound.tsv',inbound,header='bin lo hi dt v2dt')
        np.savetxt(self.output/'block_counts.tsv',np.column_stack([np.arange(BLOCKS),np.ones(BLOCKS)]),header='block count')
        np.savetxt(self.output/'orbit_class_blocks.tsv',classes,header='class block bin dt')
        (self.output/'trajectory_summary.tsv').write_text('trajectory_id rank seed block_id termination_reason t_end_s tau_in_s tau_out_s max_aphelion_km n_scatter\n'+'\n'.join(records)+'\n')

    @staticmethod
    def write_json(path: Path, value: dict) -> None:
        path.write_text(json.dumps(value))

    def run_analysis(self) -> dict:
        return analyze(self.output,read_capture_result(self.capture),make_plots=False)

    def add_population_columns(self) -> None:
        """Add full-path moments [s, km^2/s] and exhaustive incident block counts."""
        self.meta['population_bincount_version']=1
        self.meta['population_normalization']='C_geom / N_inj / shell_volume'
        self.write_json(self.output/'metadata.json',self.meta)
        data=np.loadtxt(self.output/'radial_blocks.tsv',skiprows=1)
        data[:,6]=2+np.repeat(np.arange(BLOCKS)%4,3)
        data[:,7]=data[:,6]*100
        data=np.column_stack([data,data[:,4]+data[:,10]+3,data[:,5]+data[:,11]+300])
        np.savetxt(self.output/'radial_blocks.tsv',data,header='population fixture')
        np.savetxt(self.output/'block_counts.tsv',np.column_stack([np.arange(BLOCKS),
                   np.ones(BLOCKS),np.full(BLOCKS,2),np.ones(BLOCKS)]),header='block captured injected uncaptured')

    def test_population_densities_cover_all_paths_with_common_normalization(self) -> None:
        """Check linear closure, volume-conserving rebinning and total covariance."""
        self.add_population_columns()
        result=self.run_analysis()
        data,edges,_=read_transport_blocks(self.output,self.meta)
        _,density,error=population_density(data,edges,np.full(BLOCKS,2),1e20)
        volume=4*np.pi/3*np.diff(edges**3)
        np.testing.assert_array_equal(density[:,2],density[:,:2].sum(axis=1))
        np.testing.assert_allclose(density[:,:2].T@volume,1e20/128*data[:,[12,6]].sum(axis=0))
        merged,dens,_=population_density(data,edges,np.full(BLOCKS,2),1e20,merge=3)
        self.assertIn(R_SUN_CM,merged)
        np.testing.assert_allclose(dens.T@(4*np.pi/3*np.diff(merged**3)),density.T@volume)
        total=data[:,[12,6]].sum(axis=1).reshape(BLOCKS,3)
        reps=1e20*(total.sum(axis=0)-total)/126/volume
        np.testing.assert_allclose(error[:,2],np.sqrt(63/64*((reps-reps.mean(axis=0))**2).sum(axis=0)))
        self.assertEqual(result['population_density']['N_never_captured'],64)
        annihilation=result['population_annihilation']['Gamma_s_inv']
        self.assertAlmostEqual(annihilation['total']/sum(annihilation[k] for k in ['CC','UU','CU']),1)

    def test_missing_uncaptured_counts_are_rejected(self) -> None:
        """A captured-only count file must not normalize a complete population."""
        self.add_population_columns()
        counts=np.loadtxt(self.output/'block_counts.tsv',skiprows=1)
        counts[0,3]=0
        np.savetxt(self.output/'block_counts.tsv',counts,header='block captured injected uncaptured')
        with self.assertRaisesRegex(ValueError,'cover all injections'): self.run_analysis()

    def test_old_subset_cannot_be_plotted_as_complete_population(self) -> None:
        """Legacy unscattered-only products remain usable solely for old analysis."""
        data,edges,_=read_transport_blocks(self.output,self.meta)
        with self.assertRaisesRegex(ValueError,'rerun old transport'):
            population_density(data,edges,np.full(BLOCKS,2),1e20)

    def test_all_reported_diagnostics_have_correct_delete_block_errors(self) -> None:
        result=self.run_analysis()
        self.assertEqual(set(result['central']),set(result['jackknife_se']))
        self.assertEqual(set(result['central']),set(result['jackknife_bias']))
        self.assertAlmostEqual(result['jackknife_bias_corrected']['I_tot_s2_cm3'],
                               result['central']['I_tot_s2_cm3']-result['jackknife_bias']['I_tot_s2_cm3'])
        self.assertGreater(result['jackknife_bias']['I_tot_s2_cm3'],0)
        self.assertAlmostEqual(result['central']['f_removed'],.5)
        self.assertAlmostEqual(result['jackknife_se']['f_removed'],np.sqrt(.25/63))
        reps=(672-np.array(self.cap['blocks'])[:,1])/(1280-20)
        expected=np.sqrt(63/64*((reps-reps.mean())**2).sum())
        self.assertAlmostEqual(result['jackknife_se']['C_over_C_geom'],expected)
        factor=.01*1e9/(3*K_B_EV_K*C_KM_S**2)
        reps=factor*((self.tin*self.speed2).sum()-self.tin*self.speed2)/(self.tin.sum()-self.tin)
        expected=np.sqrt(63/64*((reps-reps.mean())**2).sum())
        self.assertAlmostEqual(result['jackknife_se']['T2_in_K']/expected,1,places=10)
        summarize([self.output/'derived.json'],self.root/'scan')
        with (self.root/'scan/points.csv').open() as stream:
            row=next(csv.DictReader(stream))
        self.assertGreater(float(row['T2_in_K_se']),0)
        self.assertTrue(np.isfinite(float(row['I_tot_s2_cm3_jackknife_bias'])))

    def test_numerical_signature_mismatch_is_rejected(self) -> None:
        self.cmeta['rk_phase_tolerance']=1e-9
        self.write_json(self.capture,self.cmeta)
        with self.assertRaisesRegex(ValueError,'rk_phase_tolerance'): self.run_analysis()

    def test_rate_grid_signature_mismatch_is_rejected(self) -> None:
        # Keep the capture grid internally valid while making it incompatible.
        self.cmeta['rate_radius_points']=1000
        self.cmeta['rate_speed_points']=256
        self.cmeta['rate_max_speed']=0.02
        self.write_json(self.capture,self.cmeta)
        with self.assertRaisesRegex(ValueError,'rate_radius_points'): self.run_analysis()

    def test_capture_does_not_require_transport_radial_geometry(self) -> None:
        self.cmeta['R_remove_rsun']=550
        self.write_json(self.capture,self.cmeta)
        self.assertGreater(self.run_analysis()['central']['C_s_inv'],0)

    def test_physical_configuration_mismatch_is_rejected(self) -> None:
        self.cmeta['physical_config']={**self.cmeta['physical_config'],'DM_spin':1.0}
        self.write_json(self.capture,self.cmeta)
        with self.assertRaisesRegex(ValueError,'physical configuration'): self.run_analysis()

    def test_capture_block_mismatch_is_rejected(self) -> None:
        self.cap['blocks'][0][1]+=1
        self.write_json(self.capture,self.cap)
        with self.assertRaisesRegex(ValueError,'capture block'): self.run_analysis()

    def test_different_master_seeds_can_still_share_an_mpi_rng_stream(self) -> None:
        self.meta['mpi_ranks']=2; self.cmeta['seed']=101+1000003
        self.write_json(self.output/'metadata.json',self.meta)
        self.write_json(self.capture,self.cmeta)
        with self.assertRaisesRegex(ValueError,'disjoint'): self.run_analysis()

    def test_capture_log_and_json_have_identical_rates_and_uncertainties(self) -> None:
        log=self.root/'slurm_capture.out'
        log.write_text('scheduler header\nCAPTURE_RESULT_JSON='+json.dumps(self.cap)+'\njob finished\n')
        result=analyze(self.output,read_capture_result(log,log=True),make_plots=False)
        expected=self.run_analysis()
        for key in ['C_s_inv','Gamma_tot_s_inv','I_out_s2_cm3']:
            self.assertEqual(result['central'][key],expected['central'][key])
            self.assertEqual(result['jackknife_se'][key],expected['jackknife_se'][key])
        self.assertAlmostEqual(result['central']['C_s_inv'],1e20*672/1280)
        for option,path in [('--capture-log',log),('--capture-json',self.capture)]:
            command=[sys.executable,str(Path(__file__).with_name('analyze_point.py')),
                     str(self.output),option,str(path),'--no-plots']
            completed=subprocess.run(command,capture_output=True,text=True,timeout=15)
            self.assertEqual(completed.returncode,0,completed.stderr)

    def test_capture_log_requires_exactly_one_record(self) -> None:
        log=self.root/'capture.out'
        for text in ['no result\n',('CAPTURE_RESULT_JSON='+json.dumps(self.cap)+'\n')*2]:
            log.write_text(text)
            with self.assertRaisesRegex(ValueError,'exactly one'):
                read_capture_result(log,log=True)
        log.write_text('CAPTURE_RESULT_JSON={broken\n')
        with self.assertRaises(ValueError): read_capture_result(log,log=True)

    def test_capture_physics_and_solar_identifier_mismatches_are_rejected(self) -> None:
        for key,value in [('m_chi_GeV',.03),('sigma_SD_cm2',1e-30),('solar_model','different')]:
            with self.subTest(key=key):
                cap={**self.cap,key:value}
                with self.assertRaisesRegex(ValueError,key):
                    analyze(self.output,cap,make_plots=False)

    def test_failed_or_incomplete_capture_is_rejected(self) -> None:
        for change in [{'production_accepted':False},{'N_numerical_failures':1},
                       {'N_computational_failures':1},{'requested_samples':2000},
                       {'C_capture_s_inv':1.0},{'capture_result_schema':0}]:
            with self.subTest(change=change), self.assertRaises(ValueError):
                analyze(self.output,{**self.cap,**change},make_plots=False)

    def test_hashes_are_not_required_or_compared(self) -> None:
        self.meta.pop('source_sha256')
        self.write_json(self.output/'metadata.json',self.meta)
        self.assertGreater(self.run_analysis()['central']['C_s_inv'],0)

    def test_old_transport_schema_is_rejected(self) -> None:
        for schema in [8,9]:
            self.write_json(self.output/'metadata.json',{**self.meta,'schema_version':schema})
            with self.assertRaisesRegex(ValueError,'schema'): self.run_analysis()

    def test_halo_reference_reads_fraction_from_effective_metadata(self) -> None:
        from halo_focused_density import analyze as analyze_halo
        self.run_analysis()
        result=analyze_halo(self.output,220.,232.,544.)
        self.assertEqual(result['DM_fraction'],1.0)
        self.assertGreater(result['integral_cross_cm3'],0)

    def test_geometry_checked_in_every_block(self) -> None:
        path=self.output/'radial_blocks.tsv'; data=np.loadtxt(path,skiprows=1); data[20,2]*=1.01
        np.savetxt(path,data,header='header')
        with self.assertRaisesRegex(ValueError,'edges differ'): self.run_analysis()

    def test_block_lifetime_swap_cannot_hide_behind_global_closure(self) -> None:
        path=self.output/'radial_blocks.tsv'; data=np.loadtxt(path,skiprows=1)
        data[[0,3],4]=data[[3,0],4]
        np.savetxt(path,data,header='header')
        with self.assertRaisesRegex(ValueError,'per-block'): self.run_analysis()

    def test_incoming_path_must_close_against_its_block_histograms(self) -> None:
        path=self.output/'incident_inbound.tsv'; data=np.loadtxt(path,skiprows=1)
        data[2,3]+=1
        np.savetxt(path,data,header='header')
        with self.assertRaisesRegex(ValueError,'incident inbound block moments'): self.run_analysis()

    def test_scan_rejects_different_annihilation_coefficients_and_reused_capture(self) -> None:
        record=self.run_analysis(); other=copy.deepcopy(record)
        other['metadata']['seed']=303; other['capture_metadata']['seed']=404
        other['sigma_v_cm3_s']*=2
        path=self.root/'other.json'; self.write_json(path,other)
        with self.assertRaisesRegex(ValueError,'incompatible'):
            summarize([self.output/'derived.json',path],self.root/'scan')
        other['sigma_v_cm3_s']=record['sigma_v_cm3_s']; other['capture_metadata']['seed']=202
        self.write_json(path,other)
        with self.assertRaisesRegex(ValueError,'shared capture'):
            summarize([self.output/'derived.json',path],self.root/'scan')

    def test_missing_error_never_becomes_zero_in_scan(self) -> None:
        record=self.run_analysis(); del record['jackknife_se']['T2_in_K']
        path=self.root/'incomplete.json'; self.write_json(path,record)
        with self.assertRaisesRegex(ValueError,'uncertainty'):
            summarize([path],self.root/'scan')

    def test_accepted_marker_cannot_override_failure_counts(self) -> None:
        self.meta['N_numerical_failures']=1
        self.write_json(self.output/'metadata.json',self.meta)
        with self.assertRaisesRegex(ValueError,'failure counters'): self.run_analysis()

class GeneratorContractTests(unittest.TestCase):
    def test_generated_pairs_separate_configs_logs_and_transport_results(self) -> None:
        from prepare_transport_runs import prepare
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            manifest=prepare(root,'pilot')
            self.assertTrue((root/'manifest.json').is_file())
            self.assertTrue((root/'logs').is_dir())
            self.assertFalse((root/'results').exists())
            for run in manifest:
                text=Path(run['config']).read_text()
                self.assertIn('configs',Path(run['config']).parts)
                self.assertEqual(Path(run['config']).name,run['kind']+'.cfg')
                self.assertIn('production_mode = true;',text)
                self.assertIn('trajectory_summary_enabled = false;',text)
                self.assertIn('trajectory_events_enabled = false;',text)
                self.assertIn('trajectory_trace_rate = 0.0;',text)
                capture=run['kind']=='capture'
                self.assertIn('snapshot_enabled = '+('false' if capture else 'true')+';',text)
                self.assertEqual(run['output_dir'] is None,capture)
                if run['mass_GeV']==1.0: self.assertIn('DM_mass = 1.0;',text)

if __name__=='__main__': unittest.main()
