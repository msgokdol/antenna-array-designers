# Models and coordinate conventions

Both apps place the array or reflecting surface in the **xy plane** and use **+z** as its outward normal. Spherical **θ** is measured from +z; **φ** is measured in the xy plane from +x toward +y. Thus θ = 0° is broadside. The apps can display elevation instead, where elevation = 90° − θ for a forward-hemisphere direction.

The 2D pattern is a cut through a 3D radiation pattern, not a projection of every lobe in that pattern. A **Theta cut** varies θ at a fixed φ plane; a **Phi cut** varies φ at a fixed θ. In the reflectarray app, **Follow beam** keeps the fixed cut angle aligned with the commanded beam until unchecked. The phased-array app has an analogous cut-plane control in the Beam tab. A 3D rendering can therefore show lobes that are absent from the currently selected 2D cut.

## Phased array

The driven-array field combines element patterns with complex element excitations. Steering, amplitude taper, element rotation, and imported CST far-field data can change the result. Pattern and directivity labels depend on the selected factor and display mode. See the [phased-array walkthrough](phased-array-guide.md) for controls and CST import guidance.

## Reflectarray

For cell position **r = (x, y, 0)**, operating wavenumber **k = 2π/λ**, and target direction **u = (sinθ cosφ, sinθ sinφ)**, the synthesized reflection phase compensates incident phase and adds the outgoing phase gradient. For a horn at **F**, the incident path is **R = |r − F|** and the assigned phase is **k(R − x uₓ − y uᵧ)** modulo 2π. For a plane wave, the source angle defines an incoming direction toward the surface and the corresponding incident path is **−(x uᵢₓ + y uᵢᵧ)**.

The horn field amplitude at a cell is proportional to **max(cos α, 0)^q / R**, where α is its angle from horn boresight. Synthesized phase compensates spatial phase delay; it does not make that amplitude uniform. A normal-incidence plane wave gives uniform incident amplitude in this model. **Extra taper** adds an optional mathematical Hann weighting after illumination.

The reflected field sums the illuminated cells' complex contributions. The app normalizes that scalar field to the sum of their amplitudes. Its modeled directivity integrates power over the **upper hemisphere** and sets the field below the surface to zero. The reported dBi is therefore model directivity, **not realized gain**. Feed efficiency, spillover, physical horn pattern, polarization, cell scattering, losses, and frequency-dependent cell response are outside this model. Frequency Sweep holds cell reflection phases fixed while changing propagation frequency; it does not predict cell dispersion.

## Why matching M × N settings can look different

The phased-array default is an 8 × 8 driven array; the reflectarray default is a 20 × 20 *circularly masked* surface under an offset horn. Matching only the row and column counts still leaves shape, excitation amplitude, element pattern, phase, and normalization different. For a useful comparison, choose a full **Custom** square in both apps, match cell spacing, use **normal-incidence plane wave** with no extra taper in the reflectarray, and choose a broadside steering direction. Even then, absolute dBi labels can differ because the two apps use different radiation models.
