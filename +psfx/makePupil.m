function [A, mask, rho, theta, x, y] = makePupil(N, geom)
%MAKEPUPIL Build the generalized pupil amplitude A(x,y) on an N x N grid.
%
%   [A, mask, rho, theta, x, y] = psfx.makePupil(N, geom)
%
%   Inputs
%     N      integer grid size (N x N samples). N/q must be an integer.
%     geom   struct; missing fields get defaults:
%       .shape      'circular'|'annular'|'square'|'rect'|'hex'|'slit2'|'user'
%                   (default 'circular'; 'user' uses .userMask)
%       .obstruction epsilon in [0,0.8], inner/outer diameter ratio (annular)
%       .rectW      rect full width  as fraction of pupil diameter (default 1)
%       .rectH      rect full height as fraction of pupil diameter (default 1)
%       .slitWidth  double-slit: width of each slit, fraction of diameter (0.1)
%       .slitGap    double-slit: edge-to-edge gap, fraction of diameter (0.2)
%       .userMask   N x N amplitude map in [0,1] for shape 'user'
%       .vanes      integer 0..4 straight spider struts through the center,
%                   oriented at vaneAngle + (k-1)*180/vanes degrees (0)
%       .vaneWidth  strut width, fraction of pupil diameter (0.02)
%       .vaneAngle  orientation of the first strut, degrees (0)
%       .apod       'uniform'|'gaussian'|'cosine'|'hann'|'bessel' ('uniform')
%       .apodSigma  Gaussian sigma in units of rho/rhoMax (0.6)
%       .apodBeta   Bessel/jinc beta; 0.61 puts the first J1 zero at the edge
%       .q          padding factor N/D_pix, used only for grid spacing (4)
%
%   Outputs (units)
%     A      N x N double, amplitude transmission 0..1 (mask * apodization).
%     mask   N x N logical geometric aperture (shape + obstruction + spider),
%            BEFORE apodization. This is the support over which RMS WFE is
%            defined, so that changing apodization does not invalidate the
%            cached Zernike basis.
%     rho    N x N normalized radius; rho = 1 on the circular pupil edge.
%     theta  N x N polar angle atan2(y,x), radians.
%     x, y   N x N pupil-plane coordinates, dimensionless, pupil diameter = 2
%            and spans D_pix = N/q samples (so |x|,|y| <= 1 for a circular
%            pupil). Pixel centers: x = ((1:N) - N/2 - 1) * (2q/N), which puts
%            sample N/2+1 exactly at the pupil center (fftshift DC location).
%
%   The generalized pupil is P = A .* exp(i*2*pi*W) with W in waves; W is
%   supplied separately (see psfx.computePSF).

if nargin < 2 || isempty(geom), geom = struct(); end
if ~isfield(geom, 'q') || isempty(geom.q), geom.q = 4; end
q = geom.q;

% ---- grid coordinates (pupil diameter = 2 spans N/q samples) ----
u = ((1:N) - N/2 - 1) * (2 * q / N);
[x, y] = meshgrid(u, u);
rho = hypot(x, y);
theta = atan2(y, x);

% ---- geometric aperture ----
shape = getdef(geom, 'shape', 'circular');
switch lower(shape)
    case 'circular'
        mask = rho <= 1;
    case 'annular'
        eps0 = clampval(getdef(geom, 'obstruction', 0), 0, 0.8);
        mask = (rho <= 1) & (rho >= eps0);
    case 'square'
        mask = (abs(x) <= 1) & (abs(y) <= 1);
    case 'rect'
        w = clampval(getdef(geom, 'rectW', 1), 0.05, 1.5);
        h = clampval(getdef(geom, 'rectH', 1), 0.05, 1.5);
        mask = (abs(x) <= w) & (abs(y) <= h);
    case 'hex'
        % regular hexagon, vertices on rho = 1 (circumradius = 1)
        apothem = sqrt(3) / 2;
        mask = true(N);
        for k = 0:5
            a = deg2rad(30 + 60 * k);   % edge normals
            mask = mask & ((x * cos(a) + y * sin(a)) <= apothem);
        end
    case 'slit2'
        w = clampval(getdef(geom, 'slitWidth', 0.1), 0.01, 1);
        g = clampval(getdef(geom, 'slitGap', 0.2), 0, 2);
        ctr = (w + g) / 2;             % slit center offset from pupil center
        mask = (abs(abs(x) - ctr) <= w / 2) & (abs(y) <= 1);
    case 'user'
        um = getdef(geom, 'userMask', []);
        if isempty(um)
            error('psfx:makePupil:NoUserMask', ...
                'geom.shape = ''user'' requires geom.userMask.');
        end
        if ~isequal(size(um), [N N])
            error('psfx:makePupil:UserMaskSize', ...
                'geom.userMask must be N x N (%d x %d).', N, N);
        end
        mask = um > 0.5;
    otherwise
        error('psfx:makePupil:BadShape', 'Unknown shape ''%s''.', shape);
end


% ---- spider vanes (straight struts through the pupil center) ----
vanes = round(clampval(getdef(geom, 'vanes', 0), 0, 4));
vw = clampval(getdef(geom, 'vaneWidth', 0.02), 0, 0.25);
if vanes >= 1 && vw > 0
    va = getdef(geom, 'vaneAngle', 0);
    for k = 0:(vanes - 1)
        phi = deg2rad(va + k * 180 / vanes);
        d = abs(x * sin(phi) - y * cos(phi));   % distance to strut line
        mask = mask & ~(d <= vw);              % vw: width, pupil units (dia=2)
    end
end

if ~any(mask(:))
    error('psfx:makePupil:EmptyAperture', ...
        'Aperture definition yields an empty mask (check parameters).');
end

% ---- apodization (argument r = rho/rhoMax over the geometric mask) ----
rhoMax = max(rho(mask));
r = rho / rhoMax;
switch lower(getdef(geom, 'apod', 'uniform'))
    case 'uniform'
        apodF = 1;
    case 'gaussian'
        sg = clampval(getdef(geom, 'apodSigma', 0.6), 0.05, 5);
        apodF = exp(-(r.^2) / (2 * sg^2));
    case 'cosine'
        apodF = max(cos(pi/2 * r), 0);
    case 'hann'
        apodF = max(0.5 * (1 + cos(pi * r)), 0);
    case 'bessel'
        beta = clampval(getdef(geom, 'apodBeta', 0.61), 0.01, 5);
        den = 2 * pi * beta * r;
        apodF = ones(N);
        nz = den > 0;
        apodF(nz) = abs(2 * besselj(1, den(nz)) ./ den(nz));
    otherwise
        error('psfx:makePupil:BadApod', 'Unknown apodization.');
end

% ---- assemble amplitude ----
if strcmpi(shape, 'user')
    base = min(max(double(um), 0), 1);
else
    base = double(mask);
end
A = base .* apodF .* double(mask);

end

% ------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function v = clampval(v, lo, hi)
if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v), v = lo; end
v = min(max(v, lo), hi);
end
