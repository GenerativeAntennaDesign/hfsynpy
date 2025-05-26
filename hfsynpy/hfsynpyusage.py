from hfsynpy import Microstrip

# Example parameters (all SI units)
eps_r = 3.66  # Relative permittivity
tand = 0.0037  # Loss tangent
h = 1.524e-3  # Substrate height (1.6 mm)
t = 35e-6  # Copper thickness (35 um)
rough = 0e-6  # Surface roughness (2 um)
sigma = 1 / (1.72e-8)  # Copper conductivity (S/m)
mur = 1.0  # Relative permeability (substrate)
murc = 1.0  # Relative permeability (conductor)
frequency = 2.45e9  # Frequency (1 GHz)
z0_target = 50.0  # Target impedance (ohms)

# Create microstrip object
ms = Microstrip(
    eps_r=eps_r,
    tand=tand,
    h=h,
    t=t,
    rough=rough,
    sigma=sigma,
    mur=mur,
    murc=murc,
    frequency=frequency,
    z0_target=z0_target,
    h_top=1e20,
)

# Synthesize width for target Z0
ms.Synthesize()

# Print results
print(f"Synthesized width: {ms.width * 1e3:.4f} mm")
print(f"Effective permittivity: {ms.epsilon_eff:.4f}")
print(f"Skin depth: {ms.skin_depth * 1e6:.4f} um")
print(f"Conductor attenuation: {ms.atten_cond:.4f} dB/m")
print(f"Dielectric attenuation: {ms.atten_diel:.4f} dB/m")
