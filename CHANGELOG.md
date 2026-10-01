# CHANGELOG

## 2026-10-01 — v1.0: PupilBench PSF/OTF/MTF Explorer

- `PupilBench.m`: interactive GUI (aperture / spider / apodization / Zernike
  j=0..27 quick + table / source / display / polychromatic / image simulation /
  compare-freeze / PNG-MAT-CSV export, status bar with Strehl, RMS, FWHM,
  EE80, MTF50/10, aliasing warnings).
- `+psfx/`: `makePupil`, `zernikeTable`, `zernikeBasis`, `computePSF`,
  `computeOTF`, `encircledEnergy`, `psfMetrics`, `makeTestTarget`,
  `sumWavelengths`, `spectralWeights`.
- `tests/TestPupilBench.m`: 11 `matlab.unittest` tests, all passing.
- `README.md`, `DESIGN.md`: conventions, assumptions, screenshots list.
- Fix: `psfMetrics` edge-energy border accumulation (orientation bug).
- Fix: `uitab` uses `Title` (not `Text`); axes titles use `Title.String`;
  status labels use `.Text` (R2021a+ API correctness).
- Fix: row helpers are private methods — all call sites use `app.` prefix.
- `exportStruct` is public (Export MAT button + smoke-test diagnostics).
