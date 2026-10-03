# Screenshot gallery

These are example states of the two MATLAB apps. Numeric results depend on the array size, element pattern, illumination, phase settings, and operating frequency. Open the [phased-array guide](../phased-array-guide.md) or [reflectarray guide](../reflectarray-guide.md) to reproduce and interpret each view.

## Phased Array Designer

| View | What to look for |
| --- | --- |
| [Array and 3D overview](phased-array-overview.png) | Shape selector, wavelength-based spacing, editable element table, directivity cards, 3D pattern, and a principal 2D cut. |
| [Steered 3D pattern](phased-array-steered-3d.png) | The main lobe moves away from +z after steering. |
| [Steered 2D cut](phased-array-2d-pattern.png) | Cut-plane controls, steering φ, polar display, and absolute dBi option. |
| [Array-factor-only cut](phased-array-array-factor.png) | Pattern menu isolates the array factor from the element response. |
| [Phase map](phased-array-phase-map.png) | Per-element excitation phase and its tabular values after steering. |
| [Scan loss](phased-array-scan-loss.png) | Directivity in the commanded direction versus steering angle. |
| [Band sweep](phased-array-band-sweep.png) | Beam offset, directivity, and beamwidth versus frequency for phase shifters and true-time delay. |
| [Coverage](phased-array-coverage.png) | Scan performance across a two-dimensional region of commanded directions. |
| [Targets](phased-array-targets.png) | Optional directivity, beamwidth, sidelobe, and grating-lobe criteria. |
| [CST comparison](phased-array-cst-comparison.png) | Metric-by-metric comparison of the MATLAB cut with imported CST array data. |

### Array and 3D overview

![Phased-array geometry, controls, 3D radiation pattern, and 2D cut](phased-array-overview.png)

![Steered phased-array 3D radiation pattern](phased-array-steered-3d.png)

### Steered 2D cut and array factor

![Polar 2D pattern with steering and cut controls](phased-array-2d-pattern.png)

![Array-factor-only polar cut](phased-array-array-factor.png)

### Phase map

![Phased-array per-element phase map and table](phased-array-phase-map.png)

### Scan loss and band sweep

![Phased-array scan-loss plot](phased-array-scan-loss.png)

![Phased-array band-sweep plots and settings](phased-array-band-sweep.png)

### Coverage, targets, and CST comparison

![Phased-array scan-coverage map](phased-array-coverage.png)

![Phased-array target criteria](phased-array-targets.png)

![Phased-array comparison metrics for imported CST data](phased-array-cst-comparison.png)

## Reflectarray Designer

| View | What to look for |
| --- | --- |
| [Horn-fed geometry](reflectarray-geometry.png) | +z-facing cells, horn, design-frequency spacing, and no cell selected at launch. |
| [Plane-wave geometry](reflectarray-plane-wave.png) | One incoming red wavefront and arrow in place of the horn. |
| [Phase map](reflectarray-phase-map.png) | Incident spatial phase, progressive phase, assigned reflection phase, and continuous-aperture reference. |
| [Illumination map](reflectarray-illumination-map.png) | Incident amplitude and effective amplitude after optional taper. |
| [3D pattern](reflectarray-3d-pattern.png) | Directional reflected field over the forward hemisphere. |
| [3D pattern with cut](reflectarray-3d-cut-overlay.png) | The optional highlighted 2D cut on the 3D surface. |
| [2D pattern](reflectarray-2d-pattern.png) | Theta/phi cut, polar/Cartesian display, follow-beam control, and absolute dBi option. |
| [Frequency sweep](reflectarray-frequency-sweep.png) | Fixed-phase directivity and peak pointing as frequency changes. |

### Illumination and geometry

![Reflectarray horn-fed geometry](reflectarray-geometry.png)

![Reflectarray plane-wave geometry](reflectarray-plane-wave.png)

### Phase and amplitude

![Reflectarray four-term phase map](reflectarray-phase-map.png)

![Reflectarray illumination map](reflectarray-illumination-map.png)

### Radiation patterns

![Reflectarray 3D radiation pattern](reflectarray-3d-pattern.png)

![Reflectarray 3D pattern with selected 2D cut](reflectarray-3d-cut-overlay.png)

![Reflectarray polar 2D theta cut and steering controls](reflectarray-2d-pattern.png)

### Frequency sweep

![Reflectarray directivity and beam-pointing frequency sweep](reflectarray-frequency-sweep.png)
