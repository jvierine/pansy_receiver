import unittest
import numpy as np
from pmse.export_noise import noise_rows
from pmse.sparse_record import channel_covariance


class NoiseMigration(unittest.TestCase):
    def test_zero_detections_still_yield_all_modules_all_beams(self):
        power=np.arange(1,8)[:,None,None]*np.ones((7,5,4))
        valid=np.ones((5,4),dtype=bool); valid[0,0]=False
        power[:,0,0]=1e20  # invalid-height gate must not bias noise
        record=dict(complete=1,noise_power=power,noise_power_standard_error=power*.1,
            altitude_valid=valid,spectra=np.empty((30,7,0)),mode_start_sample=100_000_000,
            mode_stop_sample_exclusive=180_000_000)
        rows=noise_rows(record,100_000_000)
        np.testing.assert_allclose(rows['noise_power'],np.tile(np.arange(1,8),(5,1)))
        np.testing.assert_allclose(rows['times_s'],140.)

    def test_complete_marker_enforced(self):
        with self.assertRaises(ValueError): noise_rows(dict(complete=0),0)

    def test_covariance_convention_for_calibration_and_imaging(self):
        phase=np.arange(7)*.3
        z=np.tile(np.exp(1j*phase)[None,:,None],(20,1,3))
        z[::2]*=-1  # incoherent power integration, not mean voltage
        cov=channel_covariance(dict(spectra=z))
        np.testing.assert_allclose(np.angle(cov[:,1,0]),.3,atol=1e-12)
        np.testing.assert_allclose(cov,cov.conj().swapaxes(1,2),atol=1e-12)


if __name__=='__main__': unittest.main()
