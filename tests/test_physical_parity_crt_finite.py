from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def test_recorded_physical_parity_crt_finite_replay() -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts/check_physical_parity_crt_finite.py")],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert '"status": "checked"' in result.stdout
    assert '"literal_scalar_gate": "NO-GO' in result.stdout
    assert '"blockwise_global_operator_gate": "GO' in result.stdout
