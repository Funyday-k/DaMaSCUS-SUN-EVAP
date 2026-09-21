#!/usr/bin/env python3
"""Independent geometry and production-rejection tests for transport analysis."""
import json
import copy
import csv
import tempfile
import unittest
from pathlib import Path
import numpy as np
from analyze_point import analyze, BLOCKS, K_B_EV_K, C_KM_S
from summarize_transport_scan import summarize
from analyze_point import AU_CM, R_SUN_CM, chord_matrix, containment, require_accepted, angular_grid, cumulative_angle

class GeometryTests(unittest.TestCase):
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
            (p/'metadata.json').write_text(json.dumps({'schema_version':8,'workflow':'complete_captured_transport','production_accepted':False}))
            with self.assertRaises(ValueError): require_accepted(p,'complete_captured_transport')
            (p/'metadata.json').write_text(json.dumps({'schema_version':8,'workflow':'thermal_shape_validation','production_accepted':True}))
            with self.assertRaises(ValueError): require_accepted(p,'complete_captured_transport')

class AnalysisContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name); self.output=self.root/'transport'; self.capture=self.root/'capture'
        self.output.mkdir(); self.capture.mkdir()
        self.meta={'schema_version':8,'production_accepted':True,'workflow':'complete_captured_transport',
            'source_sha256':'fixture-source','seed':101,'mpi_ranks':1,'solar_model':'fixture',
            'halo_model':'SHM','halo_density_GeV_cm3':.4,'m_chi_GeV':.01,'sigma_SD_cm2':1e-32,
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
        for directory,meta in [(self.output,self.meta),(self.capture,self.cmeta)]:
            self.write_json(directory/'metadata.json',meta)
            (directory/'input.cfg').write_text('DM_mass=.01; DM_cross_section_nucleon=1e-32; DM_fraction=1; DM_distribution="SHM";')
            (directory/'solar_reference.tsv').write_text('same numerical solar reference')
        self.write_json(self.capture/'capture_summary.json',self.cap)
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
        for directory,scale in [(self.output,BLOCKS),(self.capture,1280)]:
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
        return analyze(self.output,self.capture,make_plots=False)

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
        self.cmeta['interpolation_points']=1000
        self.write_json(self.capture/'metadata.json',self.cmeta)
        with self.assertRaisesRegex(ValueError,'interpolation_points'): self.run_analysis()

    def test_rate_grid_signature_mismatch_is_rejected(self) -> None:
        # Keep the capture grid internally valid while making it incompatible.
        self.cmeta['rate_radius_points']=1000
        self.cmeta['rate_speed_points']=256
        self.cmeta['rate_max_speed']=0.02
        self.write_json(self.capture/'metadata.json',self.cmeta)
        with self.assertRaisesRegex(ValueError,'rate_radius_points'): self.run_analysis()

    def test_removal_cutoff_mismatch_is_rejected(self) -> None:
        self.cmeta['R_remove_rsun']=550
        self.write_json(self.capture/'metadata.json',self.cmeta)
        with self.assertRaisesRegex(ValueError,'R_remove_rsun'): self.run_analysis()

    def test_capture_block_mismatch_is_rejected(self) -> None:
        self.cap['blocks'][0][1]+=1
        self.write_json(self.capture/'capture_summary.json',self.cap)
        with self.assertRaisesRegex(ValueError,'capture block'): self.run_analysis()

    def test_different_master_seeds_can_still_share_an_mpi_rng_stream(self) -> None:
        self.meta['mpi_ranks']=2; self.cmeta['seed']=101+1000003
        self.write_json(self.output/'metadata.json',self.meta)
        self.write_json(self.capture/'metadata.json',self.cmeta)
        with self.assertRaisesRegex(ValueError,'disjoint'): self.run_analysis()

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

if __name__=='__main__': unittest.main()
