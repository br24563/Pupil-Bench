# DESIGN — PupilBench architecture

## State

All state in `PupilBench` private properties; no globals. `opts` is the single
source of truth (aperture, spider, apodization, 28 Zernike coeffs, source,
sampling, display, poly, sim). `cache` holds derived data, `res` the current
results, `ref` the frozen compare curves.

## Data flow

control event → `opts` (validated/clamped, inline message) → `markDirty()`
→ 50 ms coalescing `timer` → `pipeline()` → `render()` (in-place
`CData/XData/YData` updates only) → status bar + warnings + legends.

`ValueChangingFcn` on sliders sets `dirty` only; `tick()` runs at most one
`pipeline()+render()` per 50 ms (`busy` flag), so drags never queue up.

## Cached stages (`pipeline`)

| Key change | Recompute |
|---|---|
| geometry/apod (`N,q,shape,eps,rect,slit,vanes,apod,userMaskId`) | `makePupil` → `A,mask,rho,theta`; `zernikeBasis` → `B,bidx`; ideal refs (PSF/MTF/EE) |
| coefficients `wKey` | `W = B*c`, `psfStage` → PSF/OTF, metrics, sim |
| poly (`nLambda,wKind,userW`) | per-λ PSFs resampled (`sumWavelengths`) and summed |
| `lambda`/`F#` only (monochrome) | axes + metric units only — PSF samples untouched |
| display-only | `render()` mapping only (log/gamma, units, overlays) |

`maskForN` resamples a loaded user mask to `N` once per size (nearest/bilinear
via `interp2`, cached by `userMaskId`).

## Conventions

Zernike OSA/ANSI `j=(n(n+2)+m)/2` (see `+psfx/zernikeTable.m`); FFT
`fftshift(fft2(ifftshift(·)))`; PSF normalized by same-aperture ideal peak;
OTF DC = 1; PTF masked at MTF < 1e-3; pitch `λF#/q`; `df = 1/(N·pitch)`;
cutoff `1/(λF#)` — full table in `README.md`.
