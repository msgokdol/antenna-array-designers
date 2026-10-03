# Antenna Array Designers for MATLAB

Two interactive MATLAB desktop apps for exploring planar antenna arrays. **Phased Array Designer** models a driven array with steerable element excitations. **Reflectarray Designer** models a passive, phase-only reflecting surface illuminated by a horn or plane wave. Both show geometry, phase, and radiation-pattern views; they are separate apps and do not alter one another.

![Phased Array Designer overview](docs/screenshots/phased-array-overview.png)

## Start here

You need MATLAB with its desktop interface. The apps were developed and checked with **MATLAB R2025b**; other releases have not been verified. MATLAB Online and MATLAB Runtime have not been tested.

1. Download this repository using **Code → Download ZIP** on GitHub, or clone it with Git.
2. Open MATLAB and change the **Current Folder** to `phased-array` or `reflectarray` in the downloaded repository.
3. Run one of these commands in the MATLAB Command Window:

   ```matlab
   phasedArrayDesigner
   ```

   ```matlab
   reflectarrayDesigner
   ```

The phased-array folder also contains **Phased Array Designer.mlappinstall**. Double-click it in MATLAB if you prefer an app tile in the MATLAB Apps gallery. The standalone `.m` source is the clearest way to run and inspect the current code. Keep the files in each app's folder together so icons and supporting functions can be found.

## Choose an app

| | Phased Array Designer | Reflectarray Designer |
|---|---|---|
| Excitation | Driven array elements | Horn or incident plane wave reflected by a surface |
| Main controls | Array shape, element pattern, amplitude/phase, steering | Surface shape, feed/incidence, reflection phase, target beam |
| Results | 3D and 2D patterns, phase map, scan analyses, CST comparison | Geometry, four-term phase map, illumination map, 3D and 2D patterns, frequency sweep |
| Project files | MATLAB `.mat` design | MATLAB `.mat` reflectarray project |

Follow the [phased-array walkthrough](docs/phased-array-guide.md) or [reflectarray walkthrough](docs/reflectarray-guide.md) for a first experiment. The [screenshot gallery](docs/screenshots/README.md) shows the main controls and result views in both apps. The [model and coordinate notes](docs/model-and-coordinates.md) explain how to interpret angles, dB plots, and the differences between the two calculations.

![Reflectarray geometry](docs/screenshots/reflectarray-geometry.png)

## Requirements and limits

Core features of both apps use MATLAB. In Phased Array Designer, **Chebyshev** and **Taylor** tapers additionally use Signal Processing Toolbox; the app falls back to Uniform if it is unavailable. CST Studio Suite is optional and is needed only for CST data workflows.

The plots are models, not a substitute for full-wave simulation or antenna measurement. The phased-array app can use analytic or imported element patterns. The reflectarray app uses a scalar, phase-only surface model and does not predict unit-cell dispersion, reflection loss, feed spillover, coupling, polarization, or realized gain. Its optional Hann taper is a mathematical amplitude experiment, not a capability of a phase-only cell.

## Files and checks

- [`phased-array/`](phased-array/) contains the standalone MATLAB source, an optional MATLAB installer, the existing [PDF user guide](phased-array/Planar_Phased_Array_Designer_Guide.pdf), image assets, and focused checks.
- [`reflectarray/`](reflectarray/) contains the standalone MATLAB source, calculation helpers, icons, and its numerical/UI regression suite.
- [`docs/screenshots/`](docs/screenshots/README.md) contains annotated screenshots of the controls, geometry, phase and illumination maps, 2D and 3D patterns, and frequency analyses.

From each app folder in MATLAB, run `testPhasedArrayDesigner2D`, `testPhasedArrayDesignerShapes`, or `testPhasedArrayDesignerElementReference` for focused phased-array checks, and `testReflectarrayDesigner` for reflectarray numerical and UI checks. These checks require the desktop UI.

Copyright © 2026 Muhammed Said Gökdöl. **All rights reserved.** This public repository [does not grant a software reuse license](COPYRIGHT.md).
