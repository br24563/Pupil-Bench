function m = psfMetrics(PSF, otf, ctx)
%PSFMETRICS Key performance metrics of a computed PSF/OTF pair.
%
%   m = psfx.psfMetrics(PSF, otf, ctx)
%
%   Inputs
%     PSF  N x N intensity PSF, ideal peak normalized to 1
%     otf  N x N complex OTF, DC normalized to 1 (psfx.computeOTF)
%     ctx  struct with fields
%       .pitchUm    PSF sample pitch in microns (= lambda*F#/q)
%       .freqPitch  OTF frequency pitch in cycles/mm (= 1/(N*pitch), mm)
%       .cutoff     diffraction cutoff in cycles/mm (= 1/(lambda*F#))
%       .scaleUm    lambda*F# in microns (unit for the "lambda*F#" axes)
%       .lambdaNm   wavelength in nm (for RMS in nm)
%       .W          N x N wavefront error in waves (or [] to skip RMS)
%       .mask       N x N logical aperture support (RMS region)
%       .eeMode     'peak' | 'centroid' (origin for encircled energy)
%
%   Outputs (struct m, units in parentheses)
%     .strehl      peak(PSF)  (dimensionless; ideal same-aperture = 1)
%     .rmsWaves    RMS WFE over ctx.mask, piston removed (waves)
%     .rmsNm       RMS WFE * lambda (nm)
%     .marechal    exp(-(2*pi*rmsWaves)^2) Maréchal estimate (dimensionless)
%     .fwhmUm      FWHM of the PSF slice through the peak along x (microns)
%     .fwhmLf      same, in units of lambda*F# (dimensionless)
%     .ee80Um      80% encircled-energy DIAMETER (microns)
%     .ee80Lf      same, in units of lambda*F#
%     .mtf50       radial MTF = 0.5 crossing frequency (cycles/mm)
%     .mtf50Norm   same, normalized to the diffraction cutoff
%     .mtf10       radial MTF = 0.1 crossing frequency (cycles/mm)
%     .mtf10Norm   same, normalized to cutoff
%     .edgeFrac    PSF energy fraction in the 1-pixel grid border (wraparound
%                  indicator; warn when > 1e-3)
%     .nyqRatio    cutoff / Nyquist-of-PSF-sampling (>= 0.9 means the MTF
%                  reaches the sampling Nyquist limit -> warn)
%     .mtfNyqLeak  radial MTF value at the highest sampled frequency bin
%
%   Radial MTF: mean of |OTF| over 64 angles at each frequency radius;
%   crossings are found scanning outward from DC with linear interpolation.

pitchUm = ctx.pitchUm;
N = size(PSF, 1);
tot = sum(PSF(:));

% ---- Strehl & wavefront ----
m.strehl = max(PSF(:));
if ~isempty(ctx.W) && ~isempty(ctx.mask)
    wv = ctx.W(ctx.mask);
    m.rmsWaves = std(wv);
else
    m.rmsWaves = NaN;
end
m.rmsNm = m.rmsWaves * ctx.lambdaNm;
m.marechal = exp(-(2 * pi * m.rmsWaves)^2);

% ---- FWHM along x through the peak pixel ----
[~, pk] = max(PSF(:));
[pr, pc] = ind2sub([N N], pk);
s = PSF(pr, :);
half = s(pc) / 2;
L = find(s(1:pc) < half, 1, 'last');
R = find(s(pc:end) < half, 1, 'first');
fwhmPx = NaN;
if ~isempty(L) && ~isempty(R)
    R = R + pc - 1;
    dL = s(L + 1) - s(L);
    dR = s(R - 1) - s(R);
    if L >= 1 && R <= N && dL > 0 && dR > 0
        xl = L + (half - s(L)) / dL;
        xr = (R - 1) + (s(R - 1) - half) / dR;
        fwhmPx = xr - xl;
    end
end
m.fwhmUm = fwhmPx * pitchUm;
m.fwhmLf = m.fwhmUm / ctx.scaleUm;

% ---- Encircled energy 80% diameter ----
[ee, rPix] = psfx.encircledEnergy(PSF, ctx.eeMode, floor(N / 2));
i80 = find(ee >= 0.8, 1, 'first');
if isempty(i80)
    r80 = NaN;
elseif i80 == 1
    r80 = 0;
else
    r80 = rPix(i80 - 1) + (0.8 - ee(i80 - 1)) / (ee(i80) - ee(i80 - 1));
end
m.ee80Um = 2 * r80 * pitchUm;
m.ee80Lf = m.ee80Um / ctx.scaleUm;

% ---- radial MTF profile and threshold crossings ----
c = floor(N / 2) + 1;
nf = floor(N / 2) - 1;
ang = (0:63) * (2 * pi / 64);
xq = c + (0:nf)' * cos(ang);
yq = c + (0:nf)' * sin(ang);
mtf = abs(otf);
vals = interp2(mtf, xq, yq, 'linear');
prof = mean(vals, 2);
prof(1) = 1;

m.mtf50 = crossFreq(prof, 0.5) * ctx.freqPitch;
m.mtf10 = crossFreq(prof, 0.1) * ctx.freqPitch;
m.mtf50Norm = m.mtf50 / ctx.cutoff;
m.mtf10Norm = m.mtf10 / ctx.cutoff;
m.mtfNyqLeak = prof(end);

% ---- warnings ----
edgeSum = sum(PSF(1, :)) + sum(PSF(end, :)) + sum(PSF(2:end-1, 1)) ...
    + sum(PSF(2:end-1, end));
m.edgeFrac = edgeSum / tot;
m.nyqNyquist = 1000 / (2 * pitchUm);           % Nyquist, cycles/mm
m.nyqRatio = ctx.cutoff / m.nyqNyquist;

end

% ------------------------------------------------------------------
function f = crossFreq(prof, thr)
%First sub-threshold crossing scanning outward from DC, in BINS (linear
%interpolation between bins). Returns NaN if never crossed.
k = find(prof < thr, 1, 'first');
if isempty(k) || k == 1
    f = NaN;
    return
end
den = prof(k - 1) - prof(k);
if den <= 0
    f = k - 2;
else
    f = (k - 2) + (prof(k - 1) - thr) / den;
end
end
