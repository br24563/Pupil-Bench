function [ee, rPix] = encircledEnergy(PSF, mode, rMax)
%ENCIRCLEDENERGY Radial cumulative energy of the PSF about a chosen origin.
%
%   [ee, rPix] = psfx.encircledEnergy(PSF, mode, rMax)
%
%   Inputs
%     PSF    N x N intensity PSF (nonnegative; any absolute scale)
%     mode   'peak'   -> origin at the brightest pixel (default)
%            'centroid' -> intensity-weighted centroid over the whole grid
%     rMax   maximum radius in pixels (default floor(N/2)); radius samples
%            outside rMax are excluded from the accumulation, so ee(end) < 1
%            reveals energy that fell outside the inscribed circle (useful
%            for spotting wraparound/aliasing)
%   Outputs
%     ee     (rMax+1) x 1 cumulative encircled energy fraction in [0,1],
%            normalized by the total energy of the whole PSF grid
%     rPix   (rMax+1) x 1 radius vector in PIXELS, 0:rMax (1 px binning)
%
%   Units: radius in samples; convert with the PSF pixel pitch
%   pitch = lambda * F# / q  (microns) for physical axes.

if nargin < 2 || isempty(mode), mode = 'peak'; end
[Nc, Nc2] = size(PSF); %#ok<ASGLU>
N = size(PSF, 1);
if nargin < 3 || isempty(rMax), rMax = floor(N / 2); end
rMax = min(max(1, round(rMax)), floor(N * sqrt(2) / 2));

switch lower(mode)
    case 'peak'
        [~, pk] = max(PSF(:));
        [r0, c0] = ind2sub([N N], pk);
    case 'centroid'
        u = (1:N) - N/2 - 1;
        tot = sum(PSF(:));
        if tot <= 0
            error('psfx:encircledEnergy:ZeroEnergy', 'PSF has zero energy.');
        end
        c0 = 1 + sum(u .* sum(PSF, 1)) / tot;   % column (x)
        r0 = 1 + sum(u' .* sum(PSF, 2)) / tot;  % row (y)
    otherwise
        error('psfx:encircledEnergy:BadMode', 'mode must be peak|centroid.');
end

[cc, rr] = ndgrid(1:N, 1:N);
rad = round(hypot(cc - c0, rr - r0));
keep = rad <= rMax;
accum = accumarray(rad(keep) + 1, PSF(keep), [rMax + 1, 1], @sum, 0);
ee = cumsum(accum) / sum(PSF(:));
rPix = (0:rMax)';

end
