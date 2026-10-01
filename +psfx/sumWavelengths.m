function PSF = sumWavelengths(psfs, srcPitches, dstPitch, weights)
%SUMWAVELENGTHS Incoherent weighted sum of monochromatic PSFs on one grid.
%
%   PSF = psfx.sumWavelengths(psfs, srcPitches, dstPitch, weights)
%
%   Inputs
%     psfs        K x 1 cell array of N x N monochromatic intensity PSFs,
%                 each normalized so its own ideal (same aperture, W=0) peak
%                 is 1, sampled on a grid of pitch srcPitches(k)
%     srcPitches  K x 1 source pixel pitches (microns), = lambda_k * F# / q
%     dstPitch    target pixel pitch (microns); the output grid keeps size N
%                 but is relabeled/resampled to pitch dstPitch
%     weights     K x 1 spectral weights (nonnegative); normalized to sum 1
%                 inside this function
%   Outputs
%     PSF         N x N weighted sum on the destination grid (peak scale
%                 unchanged: each component's ideal peak is 1, so the ideal
%                 of the sum still peaks at 1 at the center)
%
%   Resampling uses cubic interpolation with zero outside the source grid.
%   Choosing dstPitch = min(srcPitches) guarantees every component's
%   bandlimit 1/(lambda_k*F#) stays below the destination Nyquist limit
%   q/(2*lambda_min*F#) for q >= 2, i.e. no aliasing (see README).

N = size(psfs{1}, 1);
w = weights(:) / sum(weights(:));
ax = ((1:N) - N/2 - 1);          % index offset from the DC sample
PSF = zeros(N);
for k = 1:numel(psfs)
    p = srcPitches(k);
    if abs(p - dstPitch) <= 1e-12 * dstPitch
        PSF = PSF + w(k) * psfs{k};
    else
        qIdx = ax * (dstPitch / p) + N/2 + 1;   % dest samples -> source index
        PSF = PSF + w(k) * interp2(psfs{k}, qIdx, qIdx, 'cubic', 0);
    end
end

end
