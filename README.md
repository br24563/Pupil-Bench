# PupilBench — PSF / OTF / MTF Explorer

Interactive MATLAB app: define an imaging system's pupil (shape, obstruction,
apodization, spiders, aberrations) and instantly see the PSF, OTF (MTF + PTF),
encircled energy, and Strehl ratio. Teaching tool + quick analysis tool.

Base MATLAB only (no toolboxes), R2021a or newer.

## Run

```matlab
cd('C:\Users\Brendan\OneDrive\Documents\GitHub\Pupil-Bench')
app = PupilBench;          % launch the GUI
runtests('tests')          % run the 11 unit tests (all pass, ~8 s)
```

Layout: `PupilBench.m` (single diffable class file, programmatic
`uifigure`/`uigridlayout`), `+psfx/` (pure computation functions, no UI),
`tests/TestPupilBench.m` (`matlab.unittest`), plus `DESIGN.md`, `CHANGELOG.md`.

## Conventions

| Topic | Convention |
|---|---|
| Pupil grid | `N x N`, pupil diameter spans `D_pix = N/q` pixels; `q = N/D_pix` (2, 4, 8). Pupil-plane coords: diameter = 2, `x = ((1:N)-N/2-1)*(2q/N)`, DC at sample `N/2+1` |
| Generalized pupil | `P = A .* exp(i*2*pi*W)`, `W` in waves, `A = 0` outside the aperture |
| FFT | Unnormalized DFT: `amp = fftshift(fft2(ifftshift(P)))`, `PSF = |amp|^2`. Parseval: `sum(|amp|^2) = N^2 * sum(|P|^2)` |
| PSF normalization | Divided by the ideal (same `A`, `W = 0`) peak, so the diffraction-limited unobstructed peak is exactly 1 and `max(PSF)` **is** the Strehl ratio. Note: for an annular pupil the ideal reference keeps the same obstruction, so Strehl = 1 there too; the absolute peak drop `(1-eps^2)^2` is documented in the annular test rather than folded into Strehl |
| OTF / MTF / PTF | `OTF = fftshift(fft2(ifftshift(PSF)))`, normalized to `OTF(0,0) = 1`. `MTF = |OTF|`, `PTF = angle(OTF)` set to 0 where `MTF < 1e-3` |
| Sampling | PSF pixel pitch `= lambda*F#/q`; OTF frequency pitch `= 1/(N*pitch)`; diffraction cutoff `= 1/(lambda*F#)`. Axes: microns or `lambda*F#` (PSF), cycles/mm or `f/fc` (OTF) |
| Encircled energy | Radial cumulative sum about the peak or centroid (toggle), normalized to total grid energy |
| Zernike | OSA/ANSI single index `j = (n*(n+2)+m)/2`, `j = 0..27` (n = 0..6). `m>0`: `cos(mθ)`, `m<0`: `sin(|m|θ)`. Each basis column is piston-removed and RMS-normalized, so a coefficient value in waves RMS contributes exactly that RMS for that term. Piston (j=0) is all zeros (no PSF effect). On non-circular supports modes are RMS-normalized but not orthogonal, so total RMS is measured directly (`std`), never `sqrt(sum(c^2))` |
| Units | Wavelength nm (400–1000), F# 1–20, RMS WFE in waves (and nm = waves × lambda), FWHM and EE80 diameter in µm and `λF#`, MTF50/MTF10 in cycles/mm and `f/fc` |
| Polychromatic | 3–9 wavelengths over 400–1000 nm, weights flat / D65-ish (coarse CIE D65 tabulation, linearly interpolated — an approximation) / user vector. Coefficients are waves at the λ slider; OPD is constant so `W(λ) = W·λref/λ`. Destination pitch `= λmin·F#/q` so every component's bandlimit stays below destination Nyquist for `q ≥ 2` |
| Test targets | Generated in code (USAF-style bars, Siemens star, knife edge), convolved in the frequency domain via the OTF, plus optional Gaussian noise |

## Assumptions (ambiguous / conflicting requirements)

1. **Piston j=0 kept in the table** but documented as no-effect (all-zero basis column); simpler than renumbering the OSA/ANSI table.
2. **Annular Strehl** is referenced to the same-obstruction ideal (= 1 when unaberrated); the `(1-eps²)²` absolute-peak law is asserted in the test against an unobstructed reference instead.
3. **Higher-order Zernike names** beyond the classic terms (j ≥ 15 except j=24) use generic `foilknot (n=…, x/y)` labels.
4. **D65-ish** is a 13-point coarse tabulation, not the full CIE spectrum.
5. **Metric units in poly mode** use the axis wavelength `λmin` (the destination grid); RMS nm still uses the λ slider.
6. **MATLAB `timer`** drives the 50 ms coalescing of slider drags (a documented `timer` object, not a toolbox); a `busy` flag guarantees at most one computation in flight.

## Screenshots to take

1. Diffraction-limited circular pupil: Airy PSF + analytic-matching MTF.
2. Annular (ε = 0.3) with spiders: narrowed core, higher sidelobes, MTF mid-frequency dip.
3. 0.07 waves RMS spherical: Strehl ≈ Maréchal, PTF structure visible.
4. Image-simulation tab (USAF target, blurred + noisy vs ideal).
5. Compare mode: frozen reference (dashed) vs live curve on the 1-D plots.
6. Polychromatic (5 λ, D65-ish) vs monochromatic PSF cross-section.
