function [B, idx, info] = zernikeBasis(rho, theta, mask, maxOrder)
%ZERNIKEBASIS OSA/ANSI Zernike basis evaluated on the aperture support.
%
%   [B, idx, info] = psfx.zernikeBasis(rho, theta, mask, maxOrder)
%
%   Inputs
%     rho      N x N normalized radius (from psfx.makePupil), dimensionless
%     theta    N x N polar angle, radians
%     mask     N x N logical support (the geometric aperture mask)
%     maxOrder maximum radial order n (default 6 -> terms j = 0..27)
%
%   Outputs
%     B        nPix x nJ double; column k holds term j = k-1 evaluated at the
%              support pixels, each column normalized to zero mean and unit
%              RMS (std = 1) over the support. Therefore
%                  W(support) = B * c,   c = 28x1 coefficient vector in waves
%              gives RMS(W over support) contributions of exactly |c_k| per
%              term. Piston (j = 0, zero variance) is returned as all zeros.
%     idx      nPix x 1 linear indices of the support pixels (into N x N)
%     info     struct array from psfx.zernikeTable: .j .n .m .name
%
%   Indexing: OSA/ANSI j = (n*(n+2)+m)/2 (0-based). m>0 -> cos(m*theta),
%   m<0 -> sin(|m|*theta). Over the full unit disk the modes are orthonormal;
%   over non-circular supports they are only RMS-normalized (not orthogonal),
%   so the total RMS WFE is always measured directly, never sqrt(sum(c.^2)).

if nargin < 4 || isempty(maxOrder), maxOrder = 6; end
info = psfx.zernikeTable(maxOrder);
nJ = numel(info);

idx = find(mask(:));
r = rho(idx);
th = theta(idx);
nPix = numel(idx);
B = zeros(nPix, nJ);

cacheN = -1; cacheM = -1; cacheR = [];
for k = 1:nJ
    n = info(k).n;
    mabs = abs(info(k).m);
    if n ~= cacheN || mabs ~= cacheM
        cacheR = radialPoly(n, mabs, r);
        cacheN = n; cacheM = mabs;
    end
    if info(k).m == 0
        z = cacheR;
    elseif info(k).m > 0
        z = cacheR .* cos(info(k).m * th);
    else
        z = cacheR .* sin(mabs * th);
    end
    s = std(z);
    if s < 1e-12
        B(:, k) = 0;               % piston: no effect on the PSF
    else
        B(:, k) = (z - mean(z)) / s;
    end
end

end

% ------------------------------------------------------------------
function R = radialPoly(n, mabs, rho)
%R_n^m(rho) via the standard hypergeometric sum (n <= 6, exact factorials).
R = zeros(size(rho));
kmax = floor((n - mabs) / 2);
for k = 0:kmax
    f = factorial(n - k) / (factorial(k) * ...
        factorial((n + mabs) / 2 - k) * factorial((n - mabs) / 2 - k));
    R = R + (-1)^k * f * rho.^(n - 2 * k);
end
end
