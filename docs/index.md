# hfsynpy 

A package in development for high-frequency (HF) synthesis.

This package will contain multiple classes for the synthesis of commonly used high-frequency (HF) components. Currently, only microstrip synthesis is supported.
It should produce results identical to KiCad, as the same methods are used. The implemented equations are generally applicable up to 40 GHz; however, use with caution. Always verify results with full-wave simulation tools.
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

## Attribution

This package is part of a Python translation of KiCad's C++ source code.

Original C++ code:
- © 2001 Gopal Narayanan <gopal@astro.umass.edu>
- © 2002 Claudio Girardi <claudio.girardi@ieee.org>
- © 2005, 2006 Stefan Jahn <stefan@lkcc.org>
- Modified for KiCad: 2018 Jean-Pierre Charras <jp.charras at wanadoo.fr>
- © The KiCad Developers, see AUTHORS.txt for contributors.

Python translation and modifications:
- © 2025 Dominik Mair <dominik.mair@uibk.ac.at>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.