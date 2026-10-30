import unittest
from pmse.monitor_overview import overview

class OverviewTests(unittest.TestCase):
    def test_sparse_metadata_maxima_and_gaps_are_preserved(self):
        packet=dict(axis_start_us=0,axis_stop_us=86400e6,selected_day='2026-10-30',available_days=['2026-10-30'],
            revision='test',frame_count=2,latest_mode_us=180e6,geometry={'altitude_km':[[80,81]]},
            frame_packets=[dict(start_us=0,stop_us=80e6,index=[[0,1,0],[0,2,0],[1,1,0]],snr_db=[5,9,100]),
                           dict(start_us=100e6,stop_us=180e6,index=[[0,1,1]],snr_db=[8])])
        result=overview(packet)
        self.assertEqual(result['cells'],[[0,0,9],[2,1,8]])
        self.assertEqual(result['x'],[0,80/3600,100/3600,180/3600])
        self.assertEqual(result['y'],[79.5,80.5,81.5])

if __name__=='__main__':unittest.main()
