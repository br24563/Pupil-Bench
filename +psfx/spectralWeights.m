function w = spectralWeights(kind, lambdaNm, userVec)
%SPECTRALWEIGHTS Spectral weighting samples for polychromatic PSF summation.
%
%   w = psfx.spectralWeights(kind, lambdaNm, userVec)
%
%   Inputs
%     kind      'flat' | 'd65' | 'user'
%     lambdaNm  1 x K wavelength samples in nm (the polychromatic grid)
%     userVec   1 x K user weights (used when kind = 'user'); [] otherwise
%   Outputs
%     w         1 x K raw weights >= 0 (NOT normalized; the caller normalizes,
%               see psfx.sumWavelengths)
%
%   'd65' uses a coarse CIE D65 daylight tabulation (400..1000 nm at 50 nm),
%   linearly interpolated - an approximation, labeled "D65-ish" in the UI.

lam = lambdaNm(:)';
switch lower(kind)
    case 'flat'
        w = ones(size(lam));
    case 'd65'
        t = 400:50:1000;
        v = [96 106 97 100 98 94 88 85 82 80 77 75 74] / 100;
        w = interp1(t, v, lam, 'linear', 'extrap');
        w = max(w, 0);
    case 'user'
        if isempty(userVec) || numel(userVec) ~= numel(lam)
            error('psfx:spectralWeights:UserVec', ...
                'User weight vector must have one value per wavelength.');
        end
        w = double(userVec(:))';
        w = max(w, 0);
    otherwise
        error('psfx:spectralWeights:BadKind', 'Unknown weight kind ''%s''.', kind);
end

end
