function out = computePSF(A, W, ideal)
%COMPUTEPSF Amplitude and intensity PSF of the generalized pupil.
%
%   out = psfx.computePSF(A, W, ideal)
%
%   Inputs
%     A      N x N amplitude transmission 0..1 (psfx.makePupil)
%     W      N x N wavefront error in WAVES (same grid); pass zeros(size(A))
%            or [] for the diffraction-limited case
%     ideal  optional cached ideal reference struct with fields
%              .peak  peak intensity of the ideal (W=0, same A) PSF
%              .PSF   N x N ideal intensity PSF, already normalized
%            Recomputed via FFT when omitted/empty.
%
%   Outputs (out struct)
%     .amp      N x N complex amplitude PSF, fftshift(fft2(ifftshift(P)))
%               divided by sqrt(ideal.peak) so |amp|^2 has ideal peak 1
%     .PSF      N x N intensity PSF = |amp|^2, normalized so the ideal
%               (same aperture, W = 0) peak equals 1; hence max(PSF) is the
%               Strehl ratio referenced to the same aperture
%     .strehl   max(PSF(:))
%     .ideal    the ideal reference struct used (.peak, .PSF)
%     .idealPSF same as .ideal.PSF (convenience alias)
%
%   Model: P = A .* exp(i*2*pi*W)  (W in waves); unnormalized 2-D DFT
%   (Parseval: sum(|amp_raw|^2) = N^2 * sum(|P|^2)).

if nargin < 2, W = []; end
if isempty(W), W = zeros(size(A)); elseif isscalar(W), W = W * ones(size(A)); end

if nargin < 3 || isempty(ideal)
    iAmp = fftshift(fft2(ifftshift(A)));
    ideal = struct('peak', max(abs(iAmp(:))).^2, 'PSF', abs(iAmp).^2);
    ideal.PSF = ideal.PSF / ideal.peak;
end

P = A .* exp(1i * 2 * pi * W);
amp = fftshift(fft2(ifftshift(P)));

out.amp = amp / sqrt(ideal.peak);
out.PSF = (abs(amp).^2) / ideal.peak;
out.strehl = max(out.PSF(:));
out.ideal = ideal;
out.idealPSF = ideal.PSF;

end
