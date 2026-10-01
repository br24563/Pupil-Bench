function info = zernikeTable(maxOrder)
%ZERNIKETABLE OSA/ANSI single-index Zernike term table up to radial order n.
%
%   info = psfx.zernikeTable(maxOrder)   (default maxOrder = 6)
%
%   Indexing convention: OSA/ANSI standard single index
%       j = (n*(n+2) + m)/2      (0-based; j = 0, 1, 2, ...)
%   with radial order n = 0..maxOrder and azimuthal order m = -n..n step 2.
%   For maxOrder = 6 this yields j = 0..27 (28 terms).
%
%   Inputs
%     maxOrder  maximum radial order n (0..6 supported by the basis code).
%   Outputs
%     info      1 x J struct array with fields
%       .j      OSA/ANSI index (0-based, double)
%       .n      radial order
%       .m      azimuthal order (signed)
%       .name   short descriptive name used in the UI
%
%   Sign convention of the basis functions (psfx.zernikeBasis):
%     m > 0 -> cos(m*theta),  m < 0 -> sin(|m|*theta),  m = 0 -> 1.

if nargin < 1 || isempty(maxOrder), maxOrder = 6; end
maxOrder = min(max(round(maxOrder), 0), 6);

names = { ...
    0, 'piston'; 1, 'tilt y'; 2, 'tilt x'; ...
    3, 'astig 45'; 4, 'defocus'; 5, 'astig 0'; ...
    6, 'trefoil y'; 7, 'coma y'; 8, 'coma x'; 9, 'trefoil x'; ...
   10, 'tetrafoil y'; 11, 'sec. astig 45'; 12, 'spherical (3rd)'; ...
   13, 'sec. astig 0'; 14, 'tetrafoil x'; ...
   24, 'spherical (5th)'};

info = struct('j', {}, 'n', {}, 'm', {}, 'name', {});
for n = 0:maxOrder
    for m = -n:2:n
        j = (n * (n + 2) + m) / 2;
        info(j + 1) = struct('j', j, 'n', n, 'm', m, 'name', '');
    end
end

for k = 1:numel(info)
    j = info(k).j;
    hit = find([names{:, 1}] == j, 1);
    if ~isempty(hit)
        nm = names{hit, 2};
    else
        % generic but informative label for higher-order terms
        az = abs(info(k).m);
        base = {'', 'tilt', 'astig', 'trefoil', 'tetrafoil', 'pentafoil', ...
                'hexafoil'};
        if info(k).m == 0
            nm = sprintf('spherical (n=%d)', info(k).n);
        else
            ys = 'y'; if info(k).m > 0, ys = 'x'; end
            nm = sprintf('%s (n=%d, %s)', base{az + 1}, info(k).n, ys);
        end
    end
    info(k).name = nm;
end

end
