function [img, info] = makeTestTarget(kind, N, opts)
%MAKETESTTARGET Synthetic resolution targets, generated in code (no files).
%
%   [img, info] = psfx.makeTestTarget(kind, N, opts)
%
%   Inputs
%     kind  'usaf'    USAF-style bar chart: 4 resolution levels, each with a
%                     3-bar vertical triad and a 3-bar horizontal triad
%           'siemens' Siemens star with alternating bright/dark spokes
%           'edge'    vertical knife edge (step at the image center)
%     N     grid size; use the same N as the PSF/OTF grid so the OTF can be
%           applied directly via frequency-domain convolution
%     opts  optional struct: .spokes spoke count for 'siemens' (default 36,
%           forced even so the pattern closes cleanly)
%   Outputs
%     img   N x N double in [0,1]; 1 = bright / transparent
%     info  struct: .kind, .desc (description string), .barWidths (usaf: bar
%           width in pixels per level, fine to coarse), .spokes (siemens)
%
%   All targets are computed analytically on the pixel grid.

if nargin < 3 || isempty(opts), opts = struct(); end
img = zeros(N);
info = struct('kind', kind, 'desc', '', 'barWidths', [], 'spokes', 0);

switch lower(kind)
    case 'usaf'
        % 4 rows (coarse -> fine), each: vertical triad | horizontal triad
        H = floor(N / 4);
        widths = zeros(1, 4);
        for k = 1:4
            w = max(1, round(N / 2^(k + 5)));
            widths(k) = w;
            r0 = (k - 1) * H;
            % vertical bars (left half), triad total length 5*w
            xc = round(N / 4) - floor(5 * w / 2) + 1;
            ys = r0 + max(1, round(0.15 * H)) : r0 + min(H, round(0.85 * H));
            for j = 0:2
                x1 = xc + j * 2 * w;
                img(ys, max(1, x1):min(N, x1 + w - 1)) = 1;
            end
            % horizontal bars (right half)
            x0 = round(3 * N / 4) - floor(5 * w / 2) + 1;
            xs = max(1, x0):min(N, x0 + 5 * w - 1);
            ycen = r0 + floor(H / 2) - floor(5 * w / 2) + 1;
            for j = 0:2
                y1 = ycen + j * 2 * w;
                img(max(1, y1):min(r0 + H, y1 + w - 1), xs) = 1;
            end
        end
        info.barWidths = widths;
        info.desc = sprintf(['USAF-style: 4 levels, bar widths %s px ' ...
            '(3 bars + 2 gaps per triad)'], mat2str(widths));

    case 'siemens'
        spokes = 36;
        if isfield(opts, 'spokes') && ~isempty(opts.spokes)
            spokes = opts.spokes;
        end
        spokes = max(4, 2 * round(spokes / 2));
        u = ((1:N) - (N + 1) / 2) * (2 / N);
        [X, Y] = meshgrid(u, u);
        th = atan2(Y, X);
        rad = hypot(X, Y);
        img = double((rad <= 0.95) & ...
            (mod(floor((th + pi) / (2 * pi) * spokes), 2) == 0));
        info.spokes = spokes;
        info.desc = sprintf('Siemens star, %d spokes, disk radius 0.95', spokes);

    case 'edge'
        u = ((1:N) - (N + 1) / 2) * (2 / N);
        img = double(repmat(u >= 0, N, 1));
        info.desc = 'Vertical knife edge at the image center';

    otherwise
        error('psfx:makeTestTarget:BadKind', 'Unknown target ''%s''.', kind);
end

end
