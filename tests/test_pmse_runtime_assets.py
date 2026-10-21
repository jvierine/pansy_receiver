"""Every viewer entry point must resolve one immutable JS/WASM build pair."""
import hashlib
from pathlib import Path
import re
import unittest


class RuntimeAssetsTests(unittest.TestCase):
    def test_all_imports_and_preloads_pin_matching_content_addressed_pair(self):
        root=Path(__file__).resolve().parents[1]/'pmse/rvd'
        glue=(root/'pkg/rustplotlib_web.js').read_bytes()
        wasm=(root/'pkg/rustplotlib_web_bg.wasm').read_bytes()
        release=hashlib.sha256(glue+wasm).hexdigest()[:16]
        path=root/'pkg/releases'/release
        self.assertEqual((path/'rustplotlib_web.js').read_bytes(),glue)
        self.assertEqual((path/'rustplotlib_web_bg.wasm').read_bytes(),wasm)
        for name in ('viewer.js','lib/interactive_plot.js','lib/interactive_series.js','index.html'):
            source=(root/name).read_text()
            self.assertIn(f'pkg/releases/{release}/rustplotlib_web.js',source,name)
            self.assertNotRegex(source,r'\.?\.?/pkg/rustplotlib_web\.js[\'\"]',name)
        self.assertIn("new URL('rustplotlib_web_bg.wasm', import.meta.url)",glue.decode())
        self.assertIn(f'pkg/releases/{release}/rustplotlib_web_bg.wasm',(root/'index.html').read_text())


if __name__=='__main__':unittest.main()
