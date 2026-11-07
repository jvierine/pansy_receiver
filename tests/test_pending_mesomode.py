import ast
from pathlib import Path
import unittest


source = Path(__file__).resolve().parents[1] / "quick_search_meteor.py"
node = next(n for n in ast.parse(source.read_text()).body
            if isinstance(n, ast.FunctionDef) and n.name == "pending_mesomode_blocks")
scope = {}
exec(compile(ast.Module(body=[node], type_ignores=[]), str(source), "exec"), scope)
read_pending = scope["pending_mesomode_blocks"]


class Reader:
    def __init__(self):
        self.calls = []
        self.records = {100: {"start": 10, "end": 100},
                        300: {"start": 150, "end": 300},
                        400: {"start": 350, "end": 400}}

    def read(self, lo, hi, fields):
        self.calls.append((lo, hi, fields))
        return {k: v for k, v in self.records.items() if lo <= k <= hi}


class PendingModes(unittest.TestCase):
    def test_mid_block_checkpoint_includes_spanning_block(self):
        reader = Reader()
        self.assertEqual(list(read_pending(reader, (100, 400), 200)), [300, 400])
        self.assertEqual(reader.calls, [(200, 400, ["start", "end"])])

    def test_fresh_and_caught_up(self):
        reader = Reader()
        self.assertEqual(len(read_pending(reader, (100, 400), 0)), 3)
        reader.calls.clear()
        self.assertEqual(read_pending(reader, (100, 400), 401), {})
        self.assertEqual(reader.calls, [])


if __name__ == "__main__":
    unittest.main()
