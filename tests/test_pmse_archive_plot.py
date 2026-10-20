"""Storage plots remain independent of unavailable static radar decorations."""
import copy
import io
import json
import subprocess
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from pmse.archive_plot import refresh_archive_plot,attach_archive_plot,publish_archive_plot


class ArchivePlotTests(unittest.TestCase):
    def manifest(self):
        return dict(plot_generation_status='unavailable',transfer=dict(archive=dict(
            hourly_unix_s=[0,3600],hourly_bytes=[1000,2000],
            hourly_dense_retained_fraction=[0.,.001],latest_mode_us=7_200_000_000,
            hour_tick_labels=['00:00','01:00','02:00'])))

    @staticmethod
    def render(command,**kwargs):
        assert command[1]=='--archive-only'
        out=Path(command[-1]);out.mkdir(exist_ok=True)
        for width in (640,960):(out/f'archive-growth-{width}.svg').write_text('<svg/>')

    def test_render_inventory_and_cache_unchanged_without_static_rti(self):
        with tempfile.TemporaryDirectory() as temp:
            m=self.manifest()
            with patch('pmse.archive_plot.subprocess.run',side_effect=self.render) as render:
                refresh_archive_plot(m,temp,renderer='/native-renderer')
                self.assertEqual(m['archive_plot']['status'],'ready')
                self.assertEqual(m['plot_generation_status'],'unavailable')
                original=copy.deepcopy(m['archive_plot'])
                refresh_archive_plot(m,temp,renderer='/native-renderer')
                self.assertEqual(render.call_count,1)
                self.assertEqual(m['archive_plot'],original)
                m['transfer']['archive']['hourly_bytes'][1]+=1
                refresh_archive_plot(m,temp,renderer='/native-renderer')
                self.assertEqual(render.call_count,2)
                self.assertNotEqual(m['archive_plot']['revision'],original['revision'])

    def test_failures_preserve_last_good_plots_and_do_not_block_metadata(self):
        for failure in (OSError('not installed'),subprocess.TimeoutExpired(['renderer'],120),
                        subprocess.CalledProcessError(1,['renderer'])):
            with self.subTest(failure=failure),tempfile.TemporaryDirectory() as temp:
                m=self.manifest()
                self.render(['renderer','--archive-only','spec',str(Path(temp)/'axes')])
                with patch('pmse.archive_plot.subprocess.run',side_effect=failure):
                    refresh_archive_plot(m,temp,renderer='/native-renderer')
                self.assertEqual(m['archive_plot']['status'],'delayed')
                self.assertEqual(len(m['archive_plot']['files']),2)
                self.assertIn('generated_utc',m['archive_plot'])

    def test_missing_plots_are_explicit_and_retry_after_renderer_is_available(self):
        with tempfile.TemporaryDirectory() as temp:
            m=self.manifest()
            with patch('pmse.archive_plot.subprocess.run',side_effect=FileNotFoundError):
                refresh_archive_plot(m,temp)
            self.assertEqual(m['archive_plot'],dict(status='unavailable'))
            with patch('pmse.archive_plot.subprocess.run',side_effect=self.render):
                refresh_archive_plot(m,temp)
            self.assertEqual(m['archive_plot']['status'],'ready')

    def test_web_host_attaches_analysis_status_without_running_renderer(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);m=self.manifest()
            self.render(['renderer','--archive-only','spec',str(root/'axes')])
            status=dict(status='ready',revision='analysis',files={str(w):f'axes/archive-growth-{w}.svg' for w in (640,960)})
            (root/'archive-plot-status.json').write_text(json.dumps(status))
            with patch('pmse.archive_plot.subprocess.run',side_effect=AssertionError('Periodic web plot build')):
                attach_archive_plot(m,root)
            self.assertEqual(m['archive_plot'],status)

    def test_analysis_uses_delivered_inventory_and_transfers_status_last(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            def run(command,**kwargs):
                if command[0]!='rsync':self.render(command,**kwargs)
            with patch('pmse.archive_plot.urlopen',return_value=io.BytesIO(json.dumps(self.manifest()).encode())),patch(
                    'pmse.archive_plot.subprocess.run',side_effect=run) as calls:
                publish_archive_plot('https://example/manifest.json',root,'web:data/')
            self.assertEqual(json.loads((root/'archive-plot-status.json').read_text())['status'],'ready')
            self.assertEqual(calls.call_args_list[-1].args[0][-2:],[str(root/'archive-plot-status.json'),'web:data/'])


if __name__=='__main__':unittest.main()
