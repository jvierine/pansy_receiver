"""Read-only watchdog health export. No raw reads and no receiver controls."""
from datetime import datetime
from pathlib import Path
import re


def phase_health(path):
    """Read at most 64 KiB, including rotated-log / partial-line races safely.

    Repeated 'no new pulse' does not refresh the last *measured* result's age.
    A failed/skipped measurement cannot be reported as a new healthy result.
    """
    try:
        with Path(path).open('rb') as handle:
            handle.seek(0, 2)
            size = handle.tell()
            handle.seek(max(0, size - 65536))
            lines = handle.read(65536).decode('utf-8', errors='replace').splitlines()
        for line in reversed(lines):
            stamp, _, message = line.partition(' ')
            good = message.startswith('raw phase check ok:')
            bad = ('bad raw TX phase' in message or message.startswith('raw phase check failed;'))
            unknown = (message.startswith(('raw phase reference', 'raw phase watchdog disabled', 'watchdog error:', 'raw phase check waiting'))
                       or message.startswith('raw phase check skipped:') and 'no new raw TX pulse' not in message)
            if not (good or bad or unknown):
                continue
            try:
                checked = datetime.fromisoformat(stamp).timestamp()
            except ValueError:
                continue
            # Export only the phase diagnosis, not arbitrary log command output.
            reason = message.split(': ', 1)[-1][:300]
            maximum = re.search(r'differs from reference by ([\d.]+) deg', reason)
            return dict(ok=True if good else False if bad else None, checked_unix_s=checked,
                        reason=reason, max_deviation_deg=float(maximum[1]) if maximum else None,
                        source='receiver raw-TX watchdog')
        return dict(ok=None, reason='No measured watchdog phase result in the recent log')
    except OSError:
        return dict(ok=None, reason='Watchdog phase log unavailable')
