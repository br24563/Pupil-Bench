classdef TestPupilBench < matlab.unittest.TestCase
    % Unit tests for the +psfx computation package (headless, no GUI).
    %
    % Run from the project root with:  runtests('tests')
    %
    % Covered (required) checks:
    %   1  MTF vs analytic unobstructed-circular curve, 1e-2 @ q=4, N=512
    %   2  First PSF zero at 1.22*lambda*F# (one sample) and Strehl = 1
    %   3  Annular pupil behavior + normalization choice (documented)
    %   4  0.07 waves RMS spherical: Strehl vs Marechal within 5%
    %   5  Parseval with the stated FFT-scaling convention
    %   6  MTF(0)=1, MTF real/non-negative/symmetric, PTF masked
    %   7  Sampling sanity: pitch ~ 1/q, cutoff (cyc/mm) invariant in q
    % Plus: OSA/ANSI indexing, basis RMS normalization, targets, wavelength
    % summing, spectral weights.

    properties (Constant)
        N = 512;
        Q = 4;
    end

    methods (TestClassSetup)
        function addProjectPath(~)
            here = fileparts(mfilename('fullpath'));
            addpath(fileparts(here));   % project root (parent of tests/)
        end
    end

    methods (Access = private)
        function [A, mask, rho, theta] = circularPupil(testCase, q)
            geom = struct('q', q, 'shape', 'circular');
            [A, mask, rho, theta] = psfx.makePupil(testCase.N, geom);
        end

        function r = radialMTF(~, mtf, N, q, v)
            % Radial (angle-averaged) sample of MTF at normalized
            % frequencies v = f/f_cutoff; N grid, padding q.
            c = floor(N / 2) + 1;
            rbin = v * N / q;                       % radius in freq. bins
            ang = (0:179) * pi / 90;                % 180 angles
            xq = c + rbin(:) * cos(ang);
            yq = c + rbin(:) * sin(ang);
            samp = interp2(mtf, xq, yq, 'linear');
            r = mean(samp, 2);
        end

        function r0 = firstZeroSamples(~, slice, c)
            % Radius (samples) of the first local minimum of slice(c:end)
            % (slice centered at index c) - the Airy/annular first zero.
            s = slice(c:end);
            d = diff(s);
            i = find(d(1:end-1) <= 0 & d(2:end) > 0, 1, 'first');
            if isempty(i)
                r0 = NaN;
            else
                r0 = i;      % min sits at s(i+1) -> radius = i samples
            end
        end

        function fc = measureCutoff(~, N, q, lambdaNm, fnum)
            % Last frequency (cyc/mm) where MTF along +fx stays above 1e-3.
            geom = struct('q', q);
            A = psfx.makePupil(N, geom);
            out = psfx.computePSF(A, []);
            [~, mtf] = psfx.computeOTF(out.PSF);
            c = floor(N / 2) + 1;
            pitchUm = lambdaNm * 1e-3 * fnum / q;
            freqPitch = 1000 / (N * pitchUm);       % cyc/mm per bin
            prof = mtf(c, c:end);
            k = find(prof > 1e-3, 1, 'last');
            fc = (k - 1) * freqPitch;
        end
    end

    methods (Test)
        %% ---- 1. MTF vs analytic diffraction-limited curve -------------
        function testMTFMatchesAnalytic(testCase)
            N = testCase.N; q = testCase.Q;
            A = testCase.circularPupil(q);
            out = psfx.computePSF(A, []);
            [~, mtf] = psfx.computeOTF(out.PSF);

            v = 0:0.01:0.99;
            num = testCase.radialMTF(mtf, N, q, v);
            ana = (2 / pi) * (acos(v(:)) - v(:) .* sqrt(1 - v(:).^2));
            err = max(abs(num - ana));
            testCase.verifyLessThan(err, 1e-2, ...
                sprintf('max |MTF_numeric - MTF_analytic| = %g', err));
        end

        %% ---- 2. First zero location and Strehl ------------------------
        function testFirstZeroAndStrehl(testCase)
            q = testCase.Q;
            A = testCase.circularPupil(q);
            out = psfx.computePSF(A, []);

            testCase.verifyEqual(out.strehl, 1, 'AbsTol', 1e-9);
            testCase.verifyEqual(max(out.PSF(:)), 1, 'AbsTol', 1e-9);

            c = floor(testCase.N / 2) + 1;
            r0 = testCase.firstZeroSamples(out.PSF(c, :), c);
            rExpected = 1.21966 * q;                % 1.22*lambda*F# in samples
            testCase.verifyLessThanOrEqual(abs(r0 - rExpected), 1, ...
                sprintf('first zero at %g samples, expected %g +- 1', ...
                r0, rExpected));
        end


        %% ---- 3. Annular pupil ----------------------------------------
        % Normalization choice (documented in README): the reference PSF for
        % the Strehl ratio uses the SAME aperture (spec item 4), so an
        % unaberrated annular pupil has Strehl = 1 by construction - the
        % obstruction is not penalized. The "absolute" Strehl referenced to
        % an unobstructed circular pupil of the same diameter is (1-eps^2)^2
        % (peak intensity proportional to |pupil area|^2); that is asserted
        % here end-to-end through the FFT path.
        function testAnnularBehavior(testCase)
            q = testCase.Q; N = testCase.N;
            Ac = testCase.circularPupil(q);
            outC = psfx.computePSF(Ac, []);

            geomA = struct('q', q, 'shape', 'annular', 'obstruction', 0.5);
            Aa = psfx.makePupil(N, geomA);
            outA = psfx.computePSF(Aa, []);           % annular's own reference
            outAbs = psfx.computePSF(Aa, [], outC.ideal);  % circular reference

            % (a) same-aperture normalization: unaberrated -> Strehl = 1
            testCase.verifyEqual(outA.strehl, 1, 'AbsTol', 1e-9);
            % (b) absolute Strehl vs unobstructed = (1-eps^2)^2
            eps2 = 0.5;
            testCase.verifyEqual(outAbs.strehl, (1 - eps2^2)^2, ...
                'AbsTol', 1e-2);

            c = floor(N / 2) + 1;
            % (c) narrower central lobe than the circular Airy pattern
            rCirc = testCase.firstZeroSamples(outC.PSF(c, :), c);
            rAnn = testCase.firstZeroSamples(outA.PSF(c, :), c);
            testCase.verifyLessThan(rAnn, rCirc);

            % (d) higher side lobes than Airy (1.75%)
            sideC = max(outC.PSF(c, c + rCirc + 1:end));
            sideA = max(outA.PSF(c, c + rAnn + 1:end));
            testCase.verifyGreaterThan(sideA, sideC + 0.002);
            testCase.verifyGreaterThan(sideA, 0.0175);

            % (e) mid-frequency MTF dip vs the circular pupil
            [~, oA] = psfx.computeOTF(outA.PSF);
            [~, oC] = psfx.computeOTF(outC.PSF);
            vMid = 0.35:0.05:0.60;
            mtfA = testCase.radialMTF(oA, N, q, vMid);
            mtfC = testCase.radialMTF(oC, N, q, vMid);
            testCase.verifyLessThan(min(mtfA - mtfC), -0.02, ...
                'annular MTF must sit well below circular at midband');

            % (f) encircled energy: a central obstruction pushes energy out
            % of the core into the rings - inside the core the annular EE
            % curve sits BELOW the circular one (the first zero is closer
            % in and the rings are brighter). This is the classic EE
            % penalty of a central obstruction.
            [eeA, rA] = psfx.encircledEnergy(outA.PSF, 'peak');
            [eeC, rC] = psfx.encircledEnergy(outC.PSF, 'peak');
            probe = max(2, round(0.6 * rCirc));      % inside circular first zero
            testCase.verifyLessThan(eeA(probe + 1), eeC(probe + 1));
            testCase.verifyEqual(numel(rA), numel(rC));
            testCase.verifyGreaterThan(min(eeA(end), eeC(end)), 0.98);
            % larger 80%-EE diameter: obstruction degrades energy
            % concentration (measured: 8 px vs 4 px radius at N=512, q=4)
            iA = find(eeA >= 0.8, 1, 'first');
            iC = find(eeC >= 0.8, 1, 'first');
            testCase.verifyGreaterThan(rA(iA), rC(iC));
        end

        %% ---- 4. Spherical aberration vs Marechal ----------------------
        function testSphericalMarechal(testCase)
            q = testCase.Q; N = testCase.N;
            [A, mask, rho, theta] = testCase.circularPupil(q);
            [B, idx] = psfx.zernikeBasis(rho, theta, mask);

            c = zeros(28, 1);
            c(13) = 0.07;                 % j = 12: primary spherical (3rd)
            W = zeros(N);
            W(idx) = B * c;

            % basis is unit-RMS over the aperture -> exact requested RMS
            rmsMeas = std(W(mask));
            testCase.verifyEqual(rmsMeas, 0.07, 'AbsTol', 1e-10);

            out = psfx.computePSF(A, W);
            marechal = exp(-(2 * pi * 0.07)^2);      % = 0.8241
            testCase.verifyEqual(out.strehl, marechal, 'RelTol', 0.05);
        end

        %% ---- 5. Parseval / FFT scaling convention ---------------------
        function testParseval(testCase)
            q = testCase.Q; N = testCase.N;
            [A, mask, rho, theta] = testCase.circularPupil(q);
            [B, idx] = psfx.zernikeBasis(rho, theta, mask);
            c = zeros(28, 1);
            c(5) = 0.10; c(9) = 0.05;     % defocus + coma x
            W = zeros(N); W(idx) = B * c;

            out = psfx.computePSF(A, W);
            % Convention: unnormalized DFT, PSF divided by ideal peak
            % idealPeak = |sum(A)|^2 (DC of the ideal FFT). Therefore
            %   sum(PSF) * idealPeak / N^2 = sum(|P|^2) = sum(A^2).
            lhs = sum(out.PSF(:)) * out.ideal.peak / N^2;
            rhs = sum(A(:).^2);
            testCase.verifyEqual(lhs, rhs, 'RelTol', 1e-10);
        end

        %% ---- 6. MTF properties ----------------------------------------
        function testMTFProperties(testCase)
            A = testCase.circularPupil(testCase.Q);
            out = psfx.computePSF(A, []);
            [otf, mtf, ptf] = psfx.computeOTF(out.PSF);
            N = testCase.N; c = floor(N / 2) + 1;

            testCase.verifyTrue(isreal(mtf));
            testCase.verifyGreaterThanOrEqual(min(mtf(:)), 0);
            testCase.verifyEqual(mtf(c, c), 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(imag(otf(c, c)), 0, 'AbsTol', 1e-12);

            sym = circshift(flipud(fliplr(mtf)), [1 1]);
            % centrosymmetry about (c,c): for even N the partner of index i
            % is N+2-i, hence the extra +1 shift after flipping
            testCase.verifyLessThan(max(abs(sym(:) - mtf(:))), 1e-10);

            masked = ptf(mtf < 1e-3);
            testCase.verifyTrue(all(masked == 0));
        end

        %% ---- 7. Sampling sanity ---------------------------------------
        function testSamplingSanity(testCase)
            N = testCase.N; lam = 550; fnum = 4;
            pitch2 = lam * 1e-3 * fnum / 2;
            pitch4 = lam * 1e-3 * fnum / 4;
            % spec: doubling q halves the PSF pixel pitch
            testCase.verifyEqual(pitch2 / pitch4, 2, 'RelTol', 1e-12);

            cutoff = 1000 / (lam * 1e-3 * fnum);     % cyc/mm = 1/(lambda*F#)
            fc2 = testCase.measureCutoff(N, 2, lam, fnum);
            fc4 = testCase.measureCutoff(N, 4, lam, fnum);
            % cutoff recovered from the data must match theory and be the
            % same for both padding factors
            testCase.verifyLessThan(abs(fc2 / cutoff - 1), 0.02);
            testCase.verifyLessThan(abs(fc4 / cutoff - 1), 0.02);
            testCase.verifyLessThan(abs(fc2 / fc4 - 1), 0.02);
        end

        %% ---- OSA/ANSI indexing and basis normalization ----------------
        function testZernikeIndexing(testCase)
            info = psfx.zernikeTable(6);
            testCase.verifyEqual(numel(info), 28);
            testCase.verifyEqual([info.j], 0:27);
            % spot checks of the OSA/ANSI formula j = (n(n+2)+m)/2
            testCase.verifyEqual(info(13).n, 4);   % j=12 primary spherical
            testCase.verifyEqual(info(13).m, 0);
            testCase.verifyEqual(info(5).n, 2);    % j=4 defocus
            testCase.verifyEqual(info(5).m, 0);
            testCase.verifyEqual(info(28).n, 6);   % j=27
            testCase.verifyEqual(info(28).m, 6);
            testCase.verifyEqual(info(9).m, 1);    % j=8 coma x (m=+1)

            [A, mask, rho, theta] = testCase.circularPupil(testCase.Q);
            [B, idx, infoB] = psfx.zernikeBasis(rho, theta, mask);
            testCase.verifyEqual(numel(infoB), 28);
            testCase.verifyEqual(size(B, 2), 28);
            testCase.verifyEqual(size(B, 1), nnz(mask));
            testCase.verifyEqual(numel(idx), nnz(mask));
            testCase.verifyTrue(all(B(:, 1) == 0));           % piston = 0
            for k = 2:28
                testCase.verifyEqual(std(B(:, k)), 1, 'AbsTol', 1e-9);
            end
            testCase.verifyFalse(any(isnan(B(:))));
            testCase.verifyTrue(all(A(mask) > 0));            % circular: A>0
        end

        %% ---- test target generation -----------------------------------
        function testMakeTestTarget(testCase)
            N = 128;
            kinds = {'usaf', 'siemens', 'edge'};
            for kk = 1:numel(kinds)
                [img, info] = psfx.makeTestTarget(kinds{kk}, N);
                testCase.verifyEqual(size(img), [N N]);
                testCase.verifyGreaterThanOrEqual(min(img(:)), 0);
                testCase.verifyLessThanOrEqual(max(img(:)), 1);
                testCase.verifyTrue(any(img(:) == 0) && any(img(:) == 1));
                testCase.verifyFalse(isempty(info.desc));
            end
            [~, iu] = psfx.makeTestTarget('usaf', N);
            testCase.verifyEqual(numel(iu.barWidths), 4);
            testCase.verifyTrue(all(diff(iu.barWidths) <= 0));  % coarse->fine
            edgeImg = psfx.makeTestTarget('edge', N);
            testCase.verifyTrue(all(all(edgeImg(:, 1:N/2) == 0)));
            testCase.verifyTrue(all(all(edgeImg(:, N/2+1:end) == 1)));
        end

        %% ---- polychromatic helpers ------------------------------------
        function testSumWavelengths(testCase)
            N = 64;
            p1 = rand(N); p2 = rand(N);
            w = [1, 3];
            % equal pitches -> plain normalized weighted sum (fast path)
            got = psfx.sumWavelengths({p1; p2}, [0.5; 0.5], 0.5, w);
            want = (1 * p1 + 3 * p2) / 4;
            testCase.verifyEqual(got, want, 'RelTol', 1e-12);
            % differing pitch -> resampled, same size, finite, DC-scale ok
            got2 = psfx.sumWavelengths({p1; p2}, [0.5; 1.0], 0.5, w);
            testCase.verifyEqual(size(got2), [N N]);
            testCase.verifyTrue(all(isfinite(got2(:))));
        end

        function testSpectralWeights(testCase)
            lam = 400:100:1000;
            wf = psfx.spectralWeights('flat', lam, []);
            testCase.verifyEqual(wf, ones(size(lam)));
            wd = psfx.spectralWeights('d65', lam, []);
            testCase.verifyTrue(all(wd > 0));
            testCase.verifyLessThan(max(wd) / min(wd), 2);   % smooth, gentle
            wu = psfx.spectralWeights('user', lam, ones(size(lam)));
            testCase.verifyEqual(wu, ones(size(lam)));
            testCase.verifyError(@() psfx.spectralWeights('user', lam, 1), ...
                'psfx:spectralWeights:UserVec');
        end
    end
end

