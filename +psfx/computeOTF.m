function [otf, mtf, ptf] = computeOTF(PSF)
%COMPUTEOTF Optical transfer function of a (sampled) intensity PSF.
%
%   [otf, mtf, ptf] = psfx.computeOTF(PSF)
%
%   Inputs
%     PSF   N x N intensity PSF (any nonnegative scale), DC centered as
%           produced by psfx.computePSF (peak at sample N/2+1).
%   Outputs (all N x N, fftshifted with DC at sample (N/2+1, N/2+1))
%     otf   complex OTF = fftshift(fft2(ifftshift(PSF))) / OTF(0,0),
%           normalized so otf(0) = 1 exactly
%     mtf   |otf| (dimensionless, >= 0, centrosymmetric for real PSFs)
%     ptf   angle(otf), radians in (-pi, pi], set to 0 where mtf < 1e-3
%           (phase is meaningless where the modulus vanishes)
%
%   Frequency sampling: bin k corresponds to k * df where
%   df = 1/(N * PSF_pitch) (see README conventions).

c = floor(size(PSF, 1) / 2) + 1;
otf0 = fftshift(fft2(ifftshift(PSF)));
dc = otf0(c, c);
if abs(dc) < eps
    error('psfx:computeOTF:ZeroEnergy', 'PSF has zero total energy.');
end
otf = otf0 / dc;
mtf = abs(otf);
ptf = angle(otf);
ptf(mtf < 1e-3) = 0;

end
