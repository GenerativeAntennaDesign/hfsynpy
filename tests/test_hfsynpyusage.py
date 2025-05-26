import io
import sys
import os
from hfsynpyusage import ms

def test_microstrip_output():
    captured = io.StringIO()
    sys_stdout = sys.stdout
    sys.stdout = captured
    try:
        ms.Synthesize()
        print(f"Synthesized width: {ms.width * 1e3:.3f} mm")
        print(f"Effective permittivity: {ms.epsilon_eff:.4f}")
        print(f"Skin depth: {ms.skin_depth:.3f} um")
        print(f"Conductor attenuation: {ms.atten_cond:.6f} dB/m")
        print(f"Dielectric attenuation: {ms.atten_diel:.6f} dB/m")
        print(f"Propagation delay: {ms.prop_delay:.3f} ps/cm")
    finally:
        sys.stdout = sys_stdout
    output = captured.getvalue().strip().splitlines()
    expected = [
        "Synthesized width: 3.492 mm",
        "Effective permittivity: 3.0563",
        "Skin depth: 0.660 um",
        "Conductor attenuation: 0.801981 dB/m",
        "Dielectric attenuation: 5.068878 dB/m",
        "Propagation delay: 58.315 ps/cm",
    ]
    for line, exp in zip(output, expected):
        assert line.strip() == exp, f"Expected: '{exp}', got: '{line.strip()}'"
    # Create a file to indicate test passed
    with open("test_passed.txt", "w") as f:
        f.write("All tests passed.")

if __name__ == "__main__":
    # Remove the pass file before running
    if os.path.exists("test_passed.txt"):
        os.remove("test_passed.txt")
    test_microstrip_output()
    print("All tests passed.")
