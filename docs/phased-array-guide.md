# Phased Array Designer: first experiment

![App overview](screenshots/phased-array-overview.png)

## Open the app

Use MATLAB R2025b with its desktop interface. Set MATLAB's Current Folder to `phased-array` in this repository and enter `phasedArrayDesigner`. Keep the `.png` assets beside the `.m` file. Alternatively, open **Phased Array Designer.mlappinstall** from MATLAB to add an Apps-gallery tile. The [full PDF user guide](../phased-array/Planar_Phased_Array_Designer_Guide.pdf) is included for detailed controls and CST file formats.

## Make and steer an array

1. Open the **Array** tab. Start with the default 8 × 8 uniform grid at dx = dy = 0.5 λ. M counts rows along y; N counts columns along x. Click the shape tiles or their dropdown to try Custom, Diamond, Hexagon, Octagon, Circle, and Ellipse. The center layout shows the active cells.
2. Open **Element** and choose an element pattern. Isotropic is the simplest array-factor experiment; cosine, patch, dipole, and other analytic patterns show how a single element reshapes the total beam. The CST imported-element option requires a compatible far-field file.
3. Open **Beam** and set the intended direction, for example θ = 30° and φ = 30°. θ is measured from the +z surface normal, and φ is measured from +x in the array plane. If **Auto** is off, press **Compute pattern** after changing controls.
4. Inspect the 3D radiation plot and its 2D cut. The **2D cut plane follows steering phi** option rotates a theta/elevation cut with the commanded φ. Turn it off to inspect another plane without changing the beam itself.

![Steered 2D pattern](screenshots/phased-array-2d-pattern.png)

The **2D Pattern** result offers polar and Cartesian displays and lets you plot the total pattern, array factor only, or element factor only. Choose relative dB to compare lobe shapes or **Absolute dBi** to read the model's directivity scale. A cut can miss a lobe visible in 3D if the lobe lies outside that cut plane.

For example, switching **Pattern** to **Array factor only** isolates the geometry and excitation from the element-factor shape:

![Array-factor-only 2D cut](screenshots/phased-array-array-factor.png)

## Explore phase and amplitude

The **Phase Map** result shows the phase arrangement on the aperture. The **Array** controls can change the geometry; the per-element table lets you inspect and edit amplitude, phase, and rotation. **Beam** includes amplitude taper choices. Chebyshev and Taylor taper calculations require Signal Processing Toolbox; other core operations run without it.

![Phase map](screenshots/phased-array-phase-map.png)

The **View** tab can pin a reference cut so you can compare it with a later configuration. The **Beam** tab also opens Scan Loss, Max Scan, Band Sweep, and Coverage analyses. **Scan Loss** reports directivity in the commanded direction relative to the broadside result as steering angle increases:

![Scan-loss analysis](screenshots/phased-array-scan-loss.png)

**Band Sweep** compares beam offset, commanded-direction directivity, and half-power beamwidth across frequency. At broadside the beam-offset trace is flat; steering the beam makes phase-shifter squint visible. The screen compares fixed phase shifters with true-time-delay steering:

![Band-sweep analysis](screenshots/phased-array-band-sweep.png)

Consult the PDF guide before interpreting its fixed-phase and true-time-delay choices. [Browse all app screenshots](screenshots/README.md) for the remaining main views.

**Coverage** maps scan performance across θ and φ. The **Targets** dialog lets you set optional thresholds for directivity, beamwidth, sidelobe level, and grating lobes. These are design checks, not physical tolerances.

![Scan-coverage example](screenshots/phased-array-coverage.png)

## Save, export, and compare

Use **File → Save/Save As** for a reusable MATLAB `.mat` design and **File → Open** to restore it. **File → New** resets the design. **File → Export** writes the element table as CSV or a CST-oriented array TSV. CST workflows are optional: import a supported CST element far field to use a simulated element pattern, or load a whole-array CST result for comparison with a calculated cut. The required CST angular grid and polarization conventions are described in the PDF guide.

![Example of MATLAB and imported CST cut metrics](screenshots/phased-array-cst-comparison.png)

The array and analytic element patterns are design models. They do not automatically include manufacturing tolerances, mutual coupling, feed-network loss, or mismatch. An imported CST element result may include effects contained in that imported data, but it does not make the whole array a full-wave simulation.

## Check the installation

From the `phased-array` folder, run any of the focused checks:

```matlab
testPhasedArrayDesigner2D
testPhasedArrayDesignerShapes
testPhasedArrayDesignerElementReference
```

They open the MATLAB UI and close it when done. Use `ver` to check your MATLAB release and installed products if a feature differs from the examples.
