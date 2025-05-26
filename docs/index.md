# Microstrip Transmission Line Model

---

## Overview

The `Microstrip` class models and analyzes microstrip transmission lines for PCB design, supporting both synthesis (finding trace width for a target impedance) and analysis (computing electrical properties from geometry and material parameters).

---

## Usage Example

```python
from hfsynpy import Microstrip

ms = Microstrip(
    eps_r=3.66,
    tand=0.0037,
    h=1.507e-3,
    t=35e-6,
    rough=0e-6,
    sigma=1/(1.72e-8),
    mur=1.0,
    murc=1.0,
    frequency=10e9,
    z0_target=50.0
)
ms.Synthesize()
print(f"Synthesized width: {ms.width * 1e3:.4f} mm")
print(f"Effective permittivity: {ms.epsilon_eff:.4f}")
print(f"Skin depth: {ms.skin_depth * 1e6:.4f} um")
print(f"Conductor attenuation: {ms.atten_cond:.4f} dB/m")
print(f"Dielectric attenuation: {ms.atten_diel:.4f} dB/m")
```

---

## API Reference

::: hfsynpy.Microstrip