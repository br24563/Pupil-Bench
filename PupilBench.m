classdef PupilBench < handle
    %PUPILBENCH Interactive PSF / OTF / MTF explorer (pupil -> diffraction).
    %
    %   app = PupilBench   launches the GUI (base MATLAB only, R2021a+).
    %
    %   The GUI is a thin layer: all physics lives in the +psfx package
    %   (unit-tested by tests/TestPupilBench.m). State lives entirely in
    %   object properties - there are no global variables.
    %
    %   Update path: control event -> opts (validated / clamped with an
    %   inline message) -> dirty flag -> 50 ms coalescing timer ->
    %   pipeline() (cached stages, see DESIGN.md) -> render() (updates
    %   existing graphics objects in place, never recreates them).
    %
    %   See README.md for conventions (Zernike indexing, FFT scaling,
    %   units, normalization) and DESIGN.md for the architecture.

    properties (Access = private)
        % ---- containers ----
        Fig; RootG; LeftPanel; LeftG; RightG; TabG; TabPlots; TabSim;
        % ---- control / artist / label handles ----
        ui = struct();
        ax = struct();
        h = struct();
        StatusLbl; WarnLbl; MsgLbl; TimeLbl;
        % ---- state ----
        opts;
        cache = struct();
        res = [];
        ref = [];
        zInfo = [];
        legKey = [];
        % ---- user mask ----
        userMaskRaw = [];
        userMaskN = [];
        userMaskName = '';
        userMaskId = 0;
        % ---- image-simulation target cache ----
        targetCache = struct('kind', '', 'N', -1, 'img', []);
        % ---- update coalescing ----
        dirty = false;
        busy = false;
        sync = false;
        tmr = [];
        rowCount = 0;
        lastMs = 0;
        alerted = false;
    end

    methods
        % ---- public diagnostics ----
        function S = exportStruct(app)
            %EXPORTSTRUCT Results struct (Export MAT button + smoke test).
            o = app.opts;
            r = app.res;
            cc = floor(o.N / 2) + 1;
            idx = (1:o.N) - cc;
            S = struct();
            S.created = char(datetime('now', 'Format', ...
                'yyyy-MM-dd HH:mm:ss'));
            S.opts = o;
            S.metrics = r.m;
            S.strehl = r.m.strehl;
            S.rmsWaves = r.m.rmsWaves;
            S.mtf50 = r.m.mtf50;
            S.N = o.N;
            S.PSF = r.PSF;
            S.idealPSF = r.idealPSF;
            S.OTF = r.otf;
            S.MTF = r.mtf;
            S.PTF = r.ptf;
            S.W_waves = r.W;
            S.pupilAmplitude = app.cache.A;
            S.psfAxis_um = idx * r.pitchUm;
            S.psfAxis_lambdaF = idx / o.q;
            S.freqAxis_cyc_per_mm = idx * r.freqPitch;
            S.freqAxis_norm = idx * (o.q / o.N);
            S.cutoff_cyc_per_mm = r.cutoff;
            S.encircledEnergy = r.ee;
            S.eeRadius_px = r.eeR;
            S.referenceLabel = app.refLabel();
            if ~isempty(app.ref)
                S.reference = app.ref;
            end
        end
    end

    methods
        function app = PupilBench()
            app.opts = app.defaultOpts();
            app.buildUI();
            app.updateDeps();
            t = tic;
            app.pipeline();
            app.render();
            app.lastMs = toc(t) * 1000;
            app.TimeLbl.Text = sprintf('%.1f ms', app.lastMs);
            % 50 ms coalescing timer: ValueChanging callbacks only set the
            % dirty flag, so slider drags never queue up computations.
            app.tmr = timer('ExecutionMode', 'fixedSpacing', ...
                'Period', 0.05, 'TimerFcn', @(~, ~) app.tick());
            start(app.tmr);
        end

        function delete(app)
            if ~isempty(app.tmr)
                try, stop(app.tmr); catch, end %#ok<CTCH>
                try, delete(app.tmr); catch, end %#ok<CTCH>
                app.tmr = [];
            end
            if ~isempty(app.Fig) && isvalid(app.Fig)
                delete(app.Fig);
            end
            app.Fig = [];
        end
    end

    methods (Access = private)
        % ================= options / defaults =================
        function o = defaultOpts(~)
            %DEFAULTOPTS Default parameter struct (single source of truth).
            o = struct();
            % aperture
            o.shape = 'circular';        % key matching makePupil shapes
            o.obstruction = 0;           % annular eps, 0..0.8
            o.rectW = 1; o.rectH = 1;    % rect size, fraction of diameter
            o.slitWidth = 0.1; o.slitGap = 0.2;   % double slit, fraction D
            % spider
            o.vanes = 0;                 % 0..4 straight struts
            o.vaneWidth = 0.02;          % fraction of diameter
            o.vaneAngle = 0;             % degrees
            % apodization
            o.apod = 'uniform';
            o.apodSigma = 0.6;           % Gaussian sigma, units of rhoMax
            o.apodBeta = 0.61;           % Bessel jinc beta
            % wavefront: OSA/ANSI coefficients j = 0..27, waves RMS
            o.coeff = zeros(28, 1);
            % source / sampling
            o.lambda = 550;              % nm (400..1000)
            o.fnum = 4;                  % F-number (1..20)
            o.N = 512;                   % grid size
            o.q = 4;                     % padding factor N/D_pix
            % display
            o.psfDisp = 'linear';        % linear|log|gamma (2D PSF only)
            o.gamma = 0.5;
            o.pupilView = 'amp';         % amp|phase
            o.psfUnits = 'um';           % um|lf  (microns | lambda*F#)
            o.otfUnits = 'cmm';          % cmm|norm (cyc/mm | f/fc)
            o.otfView = 'mtf';           % mtf|ptf (2D OTF axes)
            o.eeMode = 'peak';           % peak|centroid
            o.showIdeal = true;
            % polychromatic
            o.poly = false;
            o.nLambda = 5;               % 3..9 component wavelengths
            o.wKind = 'flat';            % flat|d65|user
            o.userW = [];                % user weight vector
            % image simulation
            o.target = 'usaf';           % usaf|siemens|edge
            o.noise = 0.05;              % additive gaussian sigma
        end

        % ================= layout =================
        function buildUI(app)
            app.Fig = uifigure('Name', ...
                'PupilBench - PSF / OTF / MTF Explorer', ...
                'Position', [60 50 1400 860], ...
                'CloseRequestFcn', @(~, ~) app.delete());
            app.RootG = uigridlayout(app.Fig, [1 2]);
            app.RootG.ColumnWidth = {340, '1x'};
            app.RootG.Padding = [0 0 0 0];
            app.RootG.ColumnSpacing = 0;

            % ---- left: scrollable controls ----
            app.LeftPanel = uipanel(app.RootG, 'Title', 'Controls', ...
                'FontSize', 11);
            app.place(app.LeftPanel, 1, 1);
            app.LeftG = uigridlayout(app.LeftPanel, [9 1]);
            app.LeftG.ColumnWidth = {'1x'};
            nRows = [7 3 3 8 29 4 8 4];          % content rows per section
            hSec = 30 + 26 * nRows;              % exact section heights
            app.LeftG.RowHeight = num2cell([hSec, 32]);  % + message strip
            app.LeftG.RowSpacing = 4;
            app.LeftG.Padding = [2 4 2 4];
            if isprop(app.LeftG, 'Scrollable')   % guarded: scrollable grid
                app.LeftG.Scrollable = 'on';
            end
            app.createLeft();

            % ---- right: toolbar, tabs, status ----
            app.RightG = uigridlayout(app.RootG, [4 1]);
            app.RightG.RowHeight = {40, '1x', 24, 30};
            app.RightG.Padding = [4 4 4 4];
            app.RightG.RowSpacing = 4;
            app.place(app.RightG, 1, 2);
            app.createToolbar();
            app.TabG = uitabgroup(app.RightG);
            app.place(app.TabG, 2, 1);
            app.TabPlots = uitab(app.TabG, 'Title', 'Plots');
            app.TabSim = uitab(app.TabG, 'Title', 'Image simulation');
            app.WarnLbl = uilabel(app.RightG, 'Text', '', 'FontSize', 9, ...
                'FontColor', [0.80 0.40 0], 'HorizontalAlignment', 'left', ...
                'Visible', 'off');
            app.place(app.WarnLbl, 3, 1);
            app.createStatusRow();
            app.createAxes();
            app.createSimTab();
            app.applyOptsToControls();
        end

        function g = sectionGrid(app, parent, row, title, nRows)
            %SECTIONGRID Labeled 3-column section: [label | slider | edit].
            g = uigridlayout(parent, [nRows + 1 3]);
            g.ColumnWidth = {110, '1x', 58};
            g.RowHeight = [{20}, repmat({24}, 1, nRows)];
            g.Padding = [4 4 4 6];
            g.RowSpacing = 2;
            g.ColumnSpacing = 6;
            app.place(g, row, 1);
            app.place(uilabel(g, 'Text', title, 'FontWeight', 'bold', ...
                'FontSize', 10), 1, [1 3]);
        end

        % ================= left panel =================
        function createLeft(app)
            if isempty(app.zInfo)
                app.zInfo = psfx.zernikeTable(6);
            end
            L = app.LeftG;

            % ---- 1. Aperture (7 rows) ----
            g = app.sectionGrid(L, 1, 'Aperture', 7);
            app.ui.shape = app.drow(g, 2, 'Shape', 'shape', ...
                {'Circular', 'Annular', 'Square', 'Rectangular', ...
                 'Hexagonal', 'Double slit', 'User mask'}, ...
                {'circular', 'annular', 'square', 'rect', 'hex', ...
                 'slit2', 'user'});
            app.ui.eps = app.srow(g, 3, 'Obstruction eps', ...
                'obstruction', [], 0, 0.8, '%.3f');
            app.ui.rectW = app.srow(g, 4, 'Rect width', 'rectW', ...
                [], 0.05, 1.5, '%.2f');
            app.ui.rectH = app.srow(g, 5, 'Rect height', 'rectH', ...
                [], 0.05, 1.5, '%.2f');
            app.ui.slitW = app.srow(g, 6, 'Slit width', 'slitWidth', ...
                [], 0.01, 1, '%.2f');
            app.ui.slitG = app.srow(g, 7, 'Slit gap', 'slitGap', ...
                [], 0, 1.2, '%.2f');
            app.place(uibutton(g, 'Text', 'Load mask ...', ...
                'ButtonPushedFcn', @(~, ~) app.onLoadMask()), 8, 1);
            app.ui.maskLbl = uilabel(g, 'Text', '(no user mask)', ...
                'FontSize', 9, 'HorizontalAlignment', 'left', ...
                'FontAngle', 'italic');
            app.place(app.ui.maskLbl, 8, [2 3]);

            % ---- 2. Spider / vanes (3 rows) ----
            g = app.sectionGrid(L, 2, 'Spider / vanes', 3);
            app.ui.vanes = app.drow(g, 2, 'Vanes', 'vanes', ...
                {'0', '1', '2', '3', '4'}, 0:4);
            app.ui.vaneW = app.srow(g, 3, 'Vane width', 'vaneWidth', ...
                [], 0, 0.05, '%.3f');
            app.ui.vaneA = app.srow(g, 4, 'Vane angle', 'vaneAngle', ...
                [], -180, 180, '%.0f');

            % ---- 3. Apodization (3 rows) ----
            g = app.sectionGrid(L, 3, 'Apodization', 3);
            app.ui.apod = app.drow(g, 2, 'Type', 'apod', ...
                {'Uniform', 'Gaussian', 'Cosine', 'Hann', 'Bessel (jinc)'}, ...
                {'uniform', 'gaussian', 'cosine', 'hann', 'bessel'});
            app.ui.apodS = app.srow(g, 3, 'Gaussian sigma', 'apodSigma', ...
                [], 0.1, 3, '%.2f');
            app.ui.apodB = app.srow(g, 4, 'Bessel beta', 'apodBeta', ...
                [], 0.05, 2, '%.3f');

            % ---- 4. Quick aberrations (8 rows) ----
            g = app.sectionGrid(L, 4, 'Aberrations quick access (waves RMS)', 8);
            app.ui.quickJ = [4 5 3 8 7 9 12 24];
            qlbl = {'Defocus', 'Astig 0 deg', 'Astig 45 deg', 'Coma x', ...
                'Coma y', 'Trefoil', 'Spherical 3rd', 'Spherical 5th'};
            for i = 1:8
                app.ui.quick(i) = app.srow(g, i + 1, qlbl{i}, 'coeff', ...
                    app.ui.quickJ(i), -0.5, 0.5, '%.3f');
            end

            % ---- 5. Zernike table (29 rows) ----
            g = app.sectionGrid(L, 5, ...
                'Zernike table - OSA/ANSI j = (n(n+2)+m)/2', 29);
            app.place(uilabel(g, 'Text', 'j0 piston (no effect on the PSF)', ...
                'FontSize', 8, 'FontAngle', 'italic', ...
                'HorizontalAlignment', 'left'), 2, [1 3]);
            for j = 1:27
                r = app.srow(g, j + 2, sprintf('j%d %s', j, ...
                    app.zInfo(j + 1).name), 'coeff', j, -0.5, 0.5, '%.3f');
                r.label.FontSize = 8;
                app.ui.zern(j) = r;
            end
            app.place(uibutton(g, 'Text', 'Reset all parameters', ...
                'ButtonPushedFcn', @(~, ~) app.onReset()), 30, [1 3]);

            % ---- 6. Source and sampling (4 rows) ----
            g = app.sectionGrid(L, 6, 'Source and sampling', 4);
            app.ui.lambdaRow = app.srow(g, 2, 'Wavelength (nm)', ...
                'lambda', [], 400, 1000, '%.0f');
            app.ui.fnumRow = app.srow(g, 3, 'F-number', 'fnum', ...
                [], 1, 20, '%.2f');
            app.ui.NDrop = app.drow(g, 4, 'Grid N', 'N', ...
                {'128', '256', '512', '1024'}, [128 256 512 1024]);
            app.ui.qDrop = app.drow(g, 5, 'Padding q', 'q', ...
                {'2', '4', '8'}, [2 4 8]);

            % ---- 7. Display (8 rows) ----
            g = app.sectionGrid(L, 7, 'Display', 8);
            app.ui.psfDisp = app.drow(g, 2, 'PSF display', 'psfDisp', ...
                {'Linear', 'Log', 'Gamma'}, {'linear', 'log', 'gamma'});
            app.ui.gammaRow = app.srow(g, 3, 'Gamma value', 'gamma', ...
                [], 0.1, 1, '%.2f');
            app.ui.pupilView = app.drow(g, 4, 'Pupil view', 'pupilView', ...
                {'Amplitude', 'Phase (waves)'}, {'amp', 'phase'});
            app.ui.psfUnits = app.drow(g, 5, 'PSF axis', 'psfUnits', ...
                {'microns', 'lambda*F#'}, {'um', 'lf'});
            app.ui.otfUnits = app.drow(g, 6, 'OTF axis', 'otfUnits', ...
                {'cyc/mm', 'f / cutoff'}, {'cmm', 'norm'});
            app.ui.otfView = app.drow(g, 7, '2D OTF shows', 'otfView', ...
                {'MTF', 'PTF'}, {'mtf', 'ptf'});
            app.ui.eeMode = app.drow(g, 8, 'EE origin', 'eeMode', ...
                {'Peak', 'Centroid'}, {'peak', 'centroid'});
            app.ui.idealChk = app.crow(g, 9, ...
                'Show ideal / diffraction-limit overlays', 'showIdeal');

            % ---- 8. Polychromatic (4 rows) ----
            g = app.sectionGrid(L, 8, 'Polychromatic (400-1000 nm band)', 4);
            app.ui.polyChk = app.crow(g, 2, 'Enable polychromatic sum', 'poly');
            app.ui.nLambda = app.drow(g, 3, 'Wavelengths', 'nLambda', ...
                {'3', '4', '5', '6', '7', '8', '9'}, 3:9);
            app.ui.wKind = app.drow(g, 4, 'Weights', 'wKind', ...
                {'Flat', 'D65-ish', 'User vector'}, {'flat', 'd65', 'user'});
            app.place(uilabel(g, 'Text', 'User weights', 'FontSize', 9, ...
                'HorizontalAlignment', 'left'), 5, 1);
            app.ui.userW = struct('edit', uieditfield(g, 'text', ...
                'Value', '', ...
                'ValueChangedFcn', @(~, ~) app.onUserW()));
            app.place(app.ui.userW.edit, 5, [2 3]);

            % ---- 9. Inline validation message strip ----
            mg = uigridlayout(app.LeftG, [1 1]);
            app.place(mg, 9, 1);
            mg.RowHeight = {24};
            mg.Padding = [4 4 4 4];
            app.MsgLbl = uilabel(mg, 'Text', '', 'FontSize', 9, ...
                'FontColor', [0.80 0 0], 'HorizontalAlignment', 'left', ...
                'Visible', 'off');
        end

        % ================= row helpers =================
        function r = srow(app, g, row, labelText, field, j, lo, hi, fmt)
            %SROW Slider row: [label | slider | numeric edit field].
            app.rowCount = app.rowCount + 1;
            r = struct();
            r.field = field;
            r.j = j;
            r.lims = [lo hi];
            r.fmt = fmt;
            r.id = sprintf('r%d', app.rowCount);
            r.labelText = labelText;
            r.label = uilabel(g, 'Text', labelText, 'FontSize', 9, ...
                'HorizontalAlignment', 'left');
            app.place(r.label, row, 1);
            r.slider = uislider(g, 'Limits', [lo hi], ...
                'Value', min(max(0, lo), hi));
            app.place(r.slider, row, 2);
            r.edit = uieditfield(g, 'text', 'Value', '');
            app.place(r.edit, row, 3);
            r.slider.ValueChangingFcn = @(s, e) app.onSliderEvent(r, e);
            r.slider.ValueChangedFcn = @(s, e) app.onSliderEvent(r, e);
            r.edit.ValueChangedFcn = @(s, e) app.onEditEvent(r);
        end

        function r = drow(app, g, row, labelText, field, items, itemsData)
            %DROW Dropdown row: [label | dropdown spanning 2 columns].
            app.rowCount = app.rowCount + 1;
            r = struct();
            r.field = field;
            r.id = sprintf('r%d', app.rowCount);
            r.labelText = labelText;
            r.label = uilabel(g, 'Text', labelText, 'FontSize', 9, ...
                'HorizontalAlignment', 'left');
            app.place(r.label, row, 1);
            if iscell(itemsData), v0 = itemsData{1}; else, v0 = itemsData(1); end
            r.drop = uidropdown(g, 'Items', items, 'ItemsData', itemsData, ...
                'Value', v0, ...
                'ValueChangedFcn', @(s, e) app.onDropEvent(r, e));
            app.place(r.drop, row, [2 3]);
        end

        function r = crow(app, g, row, labelText, field)
            %CROW Checkbox row spanning all three columns.
            app.rowCount = app.rowCount + 1;
            r = struct();
            r.field = field;
            r.id = sprintf('r%d', app.rowCount);
            r.labelText = labelText;
            r.chk = uicheckbox(g, 'Text', labelText, 'FontSize', 9, ...
                'ValueChangedFcn', @(s, e) app.onChkEvent(r, e));
            app.place(r.chk, row, [1 3]);
        end

        function rowEnable(~, r, on)
            %ROWENABLE Enable/disable every control of a row (inline 'off'
            % is the validation style - no dialogs).
            if on, v = 'on'; else, v = 'off'; end
            if isfield(r, 'slider'), r.slider.Enable = v; end
            if isfield(r, 'edit'),   r.edit.Enable = v; end
            if isfield(r, 'drop'),   r.drop.Enable = v; end
            if isfield(r, 'chk'),    r.chk.Enable = v; end
        end

        function setRowValue(~, r, v)
            %SETROWVALUE Push a value into both controls of a slider row.
            vv = min(max(v, r.lims(1)), r.lims(2));
            r.slider.Value = vv;
            r.edit.Value = sprintf(r.fmt, vv);
        end

        function place(~, c, row, col)
            %PLACE Position a component inside a uigridlayout.
            % row/col: scalar index, or [first last] to span a range.
            c.Layout.Row = row;
            c.Layout.Column = col;
        end

        % ================= toolbar / status =================
        function createToolbar(app)
            tb = uigridlayout(app.RightG, [1 6]);
            app.place(tb, 1, 1);
            tb.ColumnWidth = {130, 120, 10, 110, 110, 130};
            tb.Padding = [0 4 0 4];
            app.place(uibutton(tb, 'Text', 'Freeze reference', ...
                'ButtonPushedFcn', @(~, ~) app.onFreeze()), 1, 1);
            app.ui.clearBtn = uibutton(tb, 'Text', 'Clear reference', ...
                'Enable', 'off', ...
                'ButtonPushedFcn', @(~, ~) app.onClearRef());
            app.place(app.ui.clearBtn, 1, 2);
            app.place(uibutton(tb, 'Text', 'Export PNG', ...
                'ButtonPushedFcn', @(~, ~) app.onExportPNG()), 1, 4);
            app.place(uibutton(tb, 'Text', 'Export MAT', ...
                'ButtonPushedFcn', @(~, ~) app.onExportMAT()), 1, 5);
            app.place(uibutton(tb, 'Text', 'Export MTF CSV', ...
                'ButtonPushedFcn', @(~, ~) app.onExportCSV()), 1, 6);
        end

        function createStatusRow(app)
            sg = uigridlayout(app.RightG, [1 2]);
            app.place(sg, 4, 1);
            sg.ColumnWidth = {'1x', 130};
            sg.Padding = [2 2 2 2];
            app.StatusLbl = uilabel(sg, 'Text', '', 'FontSize', 9, ...
                'HorizontalAlignment', 'left');
            app.place(app.StatusLbl, 1, 1);
            app.TimeLbl = uilabel(sg, 'Text', '', 'FontSize', 9, ...
                'HorizontalAlignment', 'right', ...
                'FontColor', [0.45 0.45 0.45]);
            app.place(app.TimeLbl, 1, 2);
        end

        % ================= axes and artists =================
        function createAxes(app)
            % All artists are created ONCE; render() only updates
            % CData / XData / YData / visibility.
            g = uigridlayout(app.TabPlots, [2 3]);
            g.RowHeight = {'1x', '1x'};
            g.ColumnWidth = {'1x', '1x', '1x'};
            g.Padding = [4 4 4 4];
            g.RowSpacing = 4; g.ColumnSpacing = 4;

            % --- (1) pupil amplitude / phase ---
            ax = uiaxes(g);
            app.place(ax, 1, 1);
            ax.Box = 'on'; ax.FontSize = 9; ax.YDir = 'normal';
            axis(ax, 'image');
            ax.XLim = [-1.6 1.6]; ax.YLim = [-1.6 1.6];
            hold(ax, 'on');
            app.h.imgPupil = image(ax, [-1 1], [-1 1], zeros(2));
            colormap(ax, 'parula');
            caxis(ax, [0 1]);
            app.h.cbPupil = colorbar(ax, 'east');
            ax.XLabel.String = 'pupil units (diameter = 2)';
            ax.Title.String = 'Pupil amplitude';
            app.ax.pupil = ax;

            % --- (2) PSF image ---
            ax = uiaxes(g);
            app.place(ax, 1, 2);
            ax.Box = 'on'; ax.FontSize = 9; ax.YDir = 'normal';
            axis(ax, 'image');
            hold(ax, 'on');
            app.h.imgPSF = image(ax, [-1 1], [-1 1], zeros(2));
            colormap(ax, 'parula');
            caxis(ax, [0 1]);
            ax.XLabel.String = 'microns';
            ax.Title.String = 'PSF (linear)';
            app.ax.psf = ax;

            % --- (3) PSF cross-sections ---
            ax = uiaxes(g);
            app.place(ax, 1, 3);
            ax.Box = 'on'; ax.FontSize = 9;
            hold(ax, 'on');
            app.h.sx = line(ax, NaN, NaN, 'Color', [0 0.45 0.74]);
            app.h.sy = line(ax, NaN, NaN, 'Color', [0.85 0.33 0.10]);
            app.h.ix = line(ax, NaN, NaN, 'Color', [0 0.45 0.74], ...
                'LineStyle', '--');
            app.h.iy = line(ax, NaN, NaN, 'Color', [0.85 0.33 0.10], ...
                'LineStyle', '--');
            app.h.rx = line(ax, NaN, NaN, 'Color', [0.45 0.45 0.45], ...
                'LineStyle', ':', 'LineWidth', 1.2);
            app.h.ry = line(ax, NaN, NaN, 'Color', [0.65 0.65 0.65], ...
                'LineStyle', '-.', 'LineWidth', 1.2);
            ax.YLim = [0 1.05];
            ax.XLabel.String = 'microns';
            ax.YLabel.String = 'rel. intensity';
            ax.Title.String = 'PSF cross-sections';
            app.ax.psf1d = ax;

            app.createAxes2(g);
        end

        function createAxes2(app, g)
            % --- (4) 2D MTF/PTF with cutoff circle ---
            ax = uiaxes(g);
            app.place(ax, 2, 1);
            ax.Box = 'on'; ax.FontSize = 9; ax.YDir = 'normal';
            axis(ax, 'image');
            hold(ax, 'on');
            app.h.imgMTF = image(ax, [-1 1], [-1 1], zeros(2));
            th = linspace(0, 2 * pi, 181);
            app.h.circ = line(ax, cos(th), sin(th), 'Color', [1 0 0], ...
                'LineWidth', 1);
            colormap(ax, 'parula');
            caxis(ax, [0 1]);
            app.h.cbMtf = colorbar(ax, 'east');
            ax.XLabel.String = 'cycles/mm';
            ax.Title.String = '2D MTF (cutoff circle)';
            app.ax.mtf2d = ax;

            % --- (5) MTF cross-sections ---
            ax = uiaxes(g);
            app.place(ax, 2, 2);
            ax.Box = 'on'; ax.FontSize = 9;
            hold(ax, 'on');
            app.h.mSag = line(ax, NaN, NaN, 'Color', [0 0.45 0.74]);
            app.h.mTan = line(ax, NaN, NaN, 'Color', [0.85 0.33 0.10]);
            app.h.mDia = line(ax, NaN, NaN, 'Color', [0.93 0.69 0.13]);
            app.h.iSag = line(ax, NaN, NaN, 'Color', [0 0.45 0.74], ...
                'LineStyle', '--');
            app.h.iTan = line(ax, NaN, NaN, 'Color', [0.85 0.33 0.10], ...
                'LineStyle', '--');
            app.h.iDia = line(ax, NaN, NaN, 'Color', [0.93 0.69 0.13], ...
                'LineStyle', '--');
            app.h.rSag = line(ax, NaN, NaN, 'Color', [0.45 0.45 0.45], ...
                'LineStyle', ':', 'LineWidth', 1.2);
            app.h.rTan = line(ax, NaN, NaN, 'Color', [0.60 0.60 0.60], ...
                'LineStyle', '-.', 'LineWidth', 1.2);
            app.h.rDia = line(ax, NaN, NaN, 'Color', [0.75 0.75 0.75], ...
                'LineStyle', '--', 'LineWidth', 1.2);
            ax.YLim = [0 1.02];
            ax.XLabel.String = 'cycles/mm';
            ax.YLabel.String = 'MTF';
            ax.Title.String = 'MTF cuts: sag / tan / diag';
            app.ax.mtf1d = ax;

            % --- (6) encircled energy ---
            ax = uiaxes(g);
            app.place(ax, 2, 3);
            ax.Box = 'on'; ax.FontSize = 9;
            hold(ax, 'on');
            app.h.ee = line(ax, NaN, NaN, 'Color', [0 0.45 0.74]);
            app.h.eeI = line(ax, NaN, NaN, 'Color', [0.47 0.67 0.19], ...
                'LineStyle', '--');
            app.h.eeR = line(ax, NaN, NaN, 'Color', [0.45 0.45 0.45], ...
                'LineStyle', ':', 'LineWidth', 1.2);
            app.h.ee08 = line(ax, NaN, [0.8 0.8], 'Color', [0.5 0.5 0.5], ...
                'LineStyle', '-', 'LineWidth', 0.5);
            app.h.ee80 = line(ax, NaN, NaN, 'Color', [1 0 0], ...
                'LineStyle', ':', 'LineWidth', 1);
            ax.YLim = [0 1.02];
            ax.XLabel.String = 'radius (microns)';
            ax.YLabel.String = 'encircled energy';
            ax.Title.String = 'Encircled energy';
            app.ax.ee = ax;
        end

        function createSimTab(app)
            g = uigridlayout(app.TabSim, [2 2]);
            g.RowHeight = {58, '1x'};
            g.ColumnWidth = {'1x', '1x'};
            g.Padding = [6 6 6 6];
            g.RowSpacing = 4; g.ColumnSpacing = 6;
            gc = uigridlayout(g, [2 3]);
            app.place(gc, 1, [1 2]);
            gc.ColumnWidth = {110, '1x', 58};
            gc.RowHeight = {24, 24};
            gc.Padding = [0 2 0 2];
            gc.RowSpacing = 2; gc.ColumnSpacing = 6;
            app.place(uilabel(gc, 'Text', 'Target', 'FontSize', 9, ...
                'HorizontalAlignment', 'left'), 1, 1);
            app.ui.targetDrop = uidropdown(gc, ...
                'Items', {'USAF-style bars', 'Siemens star', 'Knife edge'}, ...
                'ItemsData', {'usaf', 'siemens', 'edge'}, 'Value', 'usaf', ...
                'ValueChangedFcn', @(~, e) app.onTargetEvent(e));
            app.place(app.ui.targetDrop, 1, [2 3]);
            app.ui.noiseRow = app.srow(gc, 2, 'Noise sigma', 'noise', ...
                [], 0, 0.3, '%.3f');

            ax = uiaxes(g);
            app.place(ax, 2, 1);
            ax.Box = 'on'; ax.FontSize = 9; ax.YDir = 'normal';
            axis(ax, 'image');
            app.h.imgTgt = image(ax, [1 10], [1 10], zeros(10));
            caxis(ax, [0 1]);
            ax.XLabel.String = 'pixels';
            ax.Title.String = 'Test target';
            app.ax.tgt = ax;

            ax = uiaxes(g);
            app.place(ax, 2, 2);
            ax.Box = 'on'; ax.FontSize = 9; ax.YDir = 'normal';
            axis(ax, 'image');
            app.h.imgSim = image(ax, [1 10], [1 10], zeros(10));
            caxis(ax, [0 1]);
            ax.XLabel.String = 'pixels';
            ax.Title.String = 'Convolved + noise';
            app.ax.sim = ax;
        end

        % ================= callbacks =================
        function onSliderEvent(app, r, e)
            % Live slider events (ValueChanging + ValueChanged): clamp to
            % the row limits, mirror into the edit field, mark dirty.
            % Actual computation happens in the coalescing timer (tick).
            v = min(max(e.Value, r.lims(1)), r.lims(2));
            r.edit.Value = sprintf(r.fmt, v);
            app.setParam(r, v);
            app.showMsg('');
        end

        function onEditEvent(app, r)
            % Numeric edit: validate; non-numeric input or out-of-range
            % values are corrected with an inline message (never a dialog).
            raw = str2double(strtrim(r.edit.Value));
            if isnan(raw)
                app.showMsg(sprintf('%s: enter a number in [%g, %g]', ...
                    r.labelText, r.lims(1), r.lims(2)));
                r.edit.Value = sprintf(r.fmt, r.slider.Value);
                return;
            end
            v = min(max(raw, r.lims(1)), r.lims(2));
            if v ~= raw
                app.showMsg(sprintf('%s clamped to %s', r.labelText, ...
                    sprintf(r.fmt, v)));
            else
                app.showMsg('');
            end
            r.edit.Value = sprintf(r.fmt, v);
            r.slider.Value = v;
            app.setParam(r, v);
        end

        function onDropEvent(app, r, e)
            v = e.Value;
            if strcmp(r.field, 'shape') && strcmp(char(v), 'user') ...
                    && isempty(app.userMaskRaw)
                app.showMsg('Load a mask file first (Load mask ...)');
                r.drop.Value = app.opts.shape;
                return;
            end
            if isnumeric(v), v = double(v); end
            app.opts.(r.field) = v;
            if ismember(r.field, {'shape', 'apod', 'poly', 'wKind'})
                app.updateDeps();
            end
            app.showMsg('');
            app.markDirty();
        end

        function onChkEvent(app, r, e)
            app.opts.(r.field) = logical(e.Value);
            if strcmp(r.field, 'poly')
                app.updateDeps();
            end
            app.markDirty();
        end

        function onTargetEvent(app, e)
            app.opts.target = char(e.Value);
            app.markDirty();
        end

        function onUserW(app)
            % User spectral weight vector: must have one nonnegative
            % value per polychromatic wavelength; otherwise the previous
            % valid vector is kept (inline message).
            v = sscanf(app.ui.userW.edit.Value, '%f')';
            if numel(v) ~= app.opts.nLambda || any(v < 0)
                app.showMsg(sprintf(['User weights: need %d ' ...
                    'nonnegative values'], app.opts.nLambda));
                app.ui.userW.edit.Value = app.weightsToString();
                return;
            end
            app.showMsg('');
            app.opts.userW = v;
            app.markDirty();
        end

        function s = weightsToString(app)
            if isempty(app.opts.userW), s = ''; return; end
            s = strtrim(sprintf('%.4g ', app.opts.userW));
        end

        function setParam(app, r, v)
            % Write a validated value into opts and mark the update dirty.
            if strcmp(r.field, 'coeff')
                app.setCoeff(r.j, v, r.id);
            else
                app.opts.(r.field) = v;
                if ismember(r.field, {'shape', 'apod', 'poly', 'wKind'})
                    app.updateDeps();
                end
            end
            app.markDirty();
        end

        function setCoeff(app, j, v, srcId)
            % Zernike coefficient: write to opts and mirror to the paired
            % controls (each coefficient has a quick-access row and a
            % table row; the source row is skipped to avoid fighting a
            % slider that the user is dragging).
            app.opts.coeff(j + 1) = v;
            if app.sync, return; end
            app.sync = true;
            zr = app.ui.zern(j);
            if ~strcmp(zr.id, srcId)
                zr.slider.Value = v;
                zr.edit.Value = sprintf(zr.fmt, v);
            end
            qi = find(app.ui.quickJ == j, 1);
            if ~isempty(qi)
                qr = app.ui.quick(qi);
                if ~strcmp(qr.id, srcId)
                    qr.slider.Value = v;
                    qr.edit.Value = sprintf(qr.fmt, v);
                end
            end
            app.sync = false;
        end

        function showMsg(app, txt)
            % Inline validation message strip (red text, no dialogs).
            if isempty(txt)
                app.MsgLbl.Visible = 'off';
            else
                app.MsgLbl.Text = txt;
                app.MsgLbl.Visible = 'on';
            end
        end

        function markDirty(app)
            app.dirty = true;   % consumed by the coalescing timer
        end

        function updateDeps(app)
            % Enable/disable controls that only apply in certain modes
            % (validation by disabling, not by error dialogs).
            o = app.opts; u = app.ui;
            app.rowEnable(u.eps,       strcmp(o.shape, 'annular'));
            app.rowEnable(u.rectW,     strcmp(o.shape, 'rect'));
            app.rowEnable(u.rectH,     strcmp(o.shape, 'rect'));
            app.rowEnable(u.slitW,     strcmp(o.shape, 'slit2'));
            app.rowEnable(u.slitG,     strcmp(o.shape, 'slit2'));
            app.rowEnable(u.apodS,     strcmp(o.apod, 'gaussian'));
            app.rowEnable(u.apodB,     strcmp(o.apod, 'bessel'));
            app.rowEnable(u.lambdaRow, ~o.poly);
            app.rowEnable(u.nLambda,   o.poly);
            app.rowEnable(u.wKind,     o.poly);
            app.rowEnable(u.userW,     o.poly && strcmp(o.wKind, 'user'));
        end

        function tick(app)
            % Coalescing timer callback: at most one computation in
            % flight; bursts of slider events collapse into one update.
            if app.busy || ~app.dirty
                return;
            end
            app.busy = true;
            app.dirty = false;
            t = tic;
            try
                app.pipeline();
                app.render();
                app.alerted = false;
            catch ME
                if ~app.alerted   % uialert only for real errors, shown once
                    app.alerted = true;
                    uialert(app.Fig, ...
                        getReport(ME, 'basic', 'hyperlinks', 'off'), ...
                        'Computation error', 'Icon', 'error');
                end
            end
            app.lastMs = toc(t) * 1000;
            app.TimeLbl.Text = sprintf('%.1f ms', app.lastMs);
            app.busy = false;
        end

        % ================= control <-> state sync =================
        function applyOptsToControls(app)
            %APPLYOPTSTOCONTROLS Push the whole opts struct into the UI
            % (initial fill and "Reset all parameters").
            o = app.opts;
            u = app.ui;
            u.shape.drop.Value = o.shape;
            app.setRowValue(u.eps, o.obstruction);
            app.setRowValue(u.rectW, o.rectW);
            app.setRowValue(u.rectH, o.rectH);
            app.setRowValue(u.slitW, o.slitWidth);
            app.setRowValue(u.slitG, o.slitGap);
            u.vanes.drop.Value = o.vanes;
            app.setRowValue(u.vaneW, o.vaneWidth);
            app.setRowValue(u.vaneA, o.vaneAngle);
            u.apod.drop.Value = o.apod;
            app.setRowValue(u.apodS, o.apodSigma);
            app.setRowValue(u.apodB, o.apodBeta);
            for j = 1:27
                app.setRowValue(u.zern(j), o.coeff(j + 1));
            end
            for i = 1:numel(u.quickJ)
                app.setRowValue(u.quick(i), o.coeff(u.quickJ(i) + 1));
            end
            app.setRowValue(u.lambdaRow, o.lambda);
            app.setRowValue(u.fnumRow, o.fnum);
            u.NDrop.drop.Value = o.N;
            u.qDrop.drop.Value = o.q;
            u.psfDisp.drop.Value = o.psfDisp;
            app.setRowValue(u.gammaRow, o.gamma);
            u.pupilView.drop.Value = o.pupilView;
            u.psfUnits.drop.Value = o.psfUnits;
            u.otfUnits.drop.Value = o.otfUnits;
            u.otfView.drop.Value = o.otfView;
            u.eeMode.drop.Value = o.eeMode;
            u.idealChk.chk.Value = o.showIdeal;
            u.polyChk.chk.Value = o.poly;
            u.nLambda.drop.Value = o.nLambda;
            u.wKind.drop.Value = o.wKind;
            u.userW.edit.Value = app.weightsToString();
            u.targetDrop.Value = o.target;
            app.setRowValue(u.noiseRow, o.noise);
        end

        function onReset(app)
            % Reset every parameter to its default (aperture, spider,
            % apodization, all Zernike terms, sampling, display).
            app.opts = app.defaultOpts();
            app.applyOptsToControls();
            app.updateDeps();
            app.showMsg('');
            app.markDirty();
        end

        function M = maskForN(app, N)
            %MASKFORN The loaded user mask resized to the current grid N.
            M = app.userMaskN;
            if ~isempty(M) && size(M, 1) == N
                return;
            end
            raw = app.userMaskRaw;
            if isempty(raw)
                M = [];
                return;
            end
            [hh, ww] = size(raw);
            [Xq, Yq] = meshgrid(linspace(1, ww, N), linspace(1, hh, N));
            M = interp2(raw, Xq, Yq, 'linear', 0);
            app.userMaskN = M;
        end

        % ================= computation pipeline =================
        function pipeline(app)
            %PIPELINE Dependency-driven, cached computation (see DESIGN.md):
            %   geom/apod change -> pupil, mask, Zernike basis, ideal ref
            %   coeff change     -> W, PSF, OTF, metrics
            %   lambda / F# only -> axes rescale + metric units (mono)
            %   poly settings    -> per-lambda PSFs resampled & summed
            o = app.opts;
            c = app.cache;
            r = app.res;
            if isempty(r), r = struct(); end

            % ---- stage keys ----
            geomKey = struct('N', o.N, 'q', o.q, 'shape', o.shape, ...
                'eps', o.obstruction, 'rw', o.rectW, 'rh', o.rectH, ...
                'sw', o.slitWidth, 'sg', o.slitGap, 'v', o.vanes, ...
                'vw', o.vaneWidth, 'va', o.vaneAngle, 'um', app.userMaskId);
            apodKey = struct('a', o.apod, 's', o.apodSigma, 'b', o.apodBeta);
            wKey = double(o.coeff(:)');
            polyKey = struct('p', o.poly, 'n', o.nLambda, 'w', o.wKind, ...
                'uw', double(o.userW(:)'));

            geomChanged = ~isfield(c, 'geomKey') || ~isequal(c.geomKey, geomKey);
            apodChanged = ~isfield(c, 'apodKey') || ~isequal(c.apodKey, apodKey);
            wChanged    = ~isfield(c, 'wKey')    || ~isequal(c.wKey, wKey);
            polyChanged = ~isfield(c, 'polyKey') || ~isequal(c.polyKey, polyKey);
            lamChanged  = ~isfield(c, 'lambdaKey') || c.lambdaKey ~= o.lambda;
            aChanged = geomChanged || apodChanged;
            idealChanged = aChanged || polyChanged || ~isfield(c, 'idealP');
            psfChanged = aChanged || wChanged || polyChanged ...
                || (o.poly && lamChanged) || ~isfield(r, 'PSF');

            % ---- 1. geometry + amplitude ----
            if aChanged
                g = struct('q', o.q, 'shape', o.shape, ...
                    'obstruction', o.obstruction, 'rectW', o.rectW, ...
                    'rectH', o.rectH, 'slitWidth', o.slitWidth, ...
                    'slitGap', o.slitGap, 'vanes', o.vanes, ...
                    'vaneWidth', o.vaneWidth, 'vaneAngle', o.vaneAngle, ...
                    'apod', o.apod, 'apodSigma', o.apodSigma, ...
                    'apodBeta', o.apodBeta, 'userMask', app.maskForN(o.N));
                try
                    [A, mask, rho, theta] = psfx.makePupil(o.N, g);
                catch ME
                    % input-level problem (e.g. empty user mask): fall
                    % back inline with a message, not a dialog
                    app.showMsg(['Aperture: ' ME.message ...
                        ' - using circular.']);
                    app.opts.shape = 'circular';
                    app.ui.shape.drop.Value = 'circular';
                    g.shape = 'circular';
                    g.userMask = [];
                    [A, mask, rho, theta] = psfx.makePupil(o.N, g);
                end
                c.A = A; c.mask = mask; c.rho = rho; c.theta = theta;
            end

            % ---- 2. Zernike basis (depends on mask only, not apod) ----
            if geomChanged
                [c.B, c.bidx] = psfx.zernikeBasis(c.rho, c.theta, c.mask);
            end

            % ---- 3. ideal reference (W = 0, same A; poly-resampled) ----
            if idealChanged
                tmp = psfx.computePSF(c.A, []);
                c.idealP = tmp.ideal;
                if o.poly
                    [~, wts, pitches, dstPitch] = app.polyGrid();
                    K = numel(pitches);
                    c.idealDisp = psfx.sumWavelengths( ...
                        repmat({c.idealP.PSF}, K, 1), pitches, dstPitch, wts);
                else
                    c.idealDisp = c.idealP.PSF;
                end
                c.idealDisp = c.idealDisp / max(c.idealDisp(:));
                [c.idealEE, c.idealEER] = ...
                    psfx.encircledEnergy(c.idealDisp, 'peak');
                [~, c.idealMtf] = psfx.computeOTF(c.idealDisp);
            end

            % ---- store keys for the next incremental run ----
            c.geomKey = geomKey;
            c.apodKey = apodKey;
            c.wKey = wKey;
            c.polyKey = polyKey;
            c.lambdaKey = o.lambda;
            app.res = r;
            app.cache = c;
            if psfChanged
                app.psfStage();
            end
            app.metricsStage();
            app.simStage();
        end

        function [lam, wts, pitches, dstPitch] = polyGrid(app)
            %POLYGRID Polychromatic sampling grid (README conventions):
            % band 400-1000 nm, destination pitch = lambda_min*F#/q, so
            % every component's bandlimit stays below destination Nyquist.
            o = app.opts;
            lam = linspace(400, 1000, o.nLambda);
            wts = psfx.spectralWeights(o.wKind, lam, o.userW);
            wts = wts / sum(wts);
            pitches = lam * 1e-3 * o.fnum / o.q;
            dstPitch = 0.4e-3 * o.fnum / o.q;
        end

        function psfStage(app)
            %PSFSTAGE W -> (mono or per-lambda) PSF -> OTF/MTF/PTF.
            o = app.opts;
            c = app.cache;
            r = app.res;

            W = zeros(o.N);
            if any(o.coeff ~= 0)
                W(c.bidx) = c.B * o.coeff;
            end

            if ~o.poly
                out = psfx.computePSF(c.A, W, c.idealP);
                psf = out.PSF;
            else
                % coefficients are waves at the lambda slider (lambdaRef);
                % optical OPD is constant -> W(lambda) = W * lamRef/lam
                [lam, wts, pitches, dstPitch] = app.polyGrid();
                cells = cell(numel(lam), 1);
                for k = 1:numel(lam)
                    ok = psfx.computePSF(c.A, W * (o.lambda / lam(k)), ...
                        c.idealP);
                    cells{k} = ok.PSF;
                end
                psf = psfx.sumWavelengths(cells, pitches, dstPitch, wts);
            end

            ipk = max(c.idealDisp(:));
            r.PSF = psf / ipk;
            r.idealPSF = c.idealDisp / ipk;    % ideal peak == 1 exactly
            [r.otf, r.mtf, r.ptf] = psfx.computeOTF(r.PSF);
            r.W = W;
            app.res = r;
        end

        function metricsStage(app)
            %METRICSSTAGE Strehl, RMS WFE, FWHM, EE80, MTF50/10, warnings.
            % Axis scales depend on lambda/F#/q but the PSF samples do
            % not (mono), so this stage runs on every update cheaply.
            o = app.opts;
            r = app.res;
            c = app.cache;
            if o.poly, lamAx = 0.4e-3; else, lamAx = o.lambda * 1e-3; end
            pitchUm = lamAx * o.fnum / o.q;         % PSF pixel pitch, um
            ctx = struct('pitchUm', pitchUm, ...
                'freqPitch', 1000 / (o.N * pitchUm), ...   % cyc/mm
                'cutoff', 1000 / (lamAx * o.fnum), ...     % cyc/mm
                'scaleUm', lamAx * o.fnum, ...             % lambda*F#, um
                'lambdaNm', o.lambda, 'W', r.W, 'mask', c.mask, ...
                'eeMode', o.eeMode);
            r.m = psfx.psfMetrics(r.PSF, r.otf, ctx);
            [r.ee, r.eeR] = psfx.encircledEnergy(r.PSF, o.eeMode, ...
                floor(o.N / 2));
            r.pitchUm = pitchUm;
            r.freqPitch = ctx.freqPitch;
            r.cutoff = ctx.cutoff;
            r.scaleUm = ctx.scaleUm;
            app.res = r;
        end

        function simStage(app)
            %SIMSTAGE Convolve the test target with the PSF in the
            % frequency domain via the OTF, then add optional noise.
            o = app.opts;
            r = app.res;
            tg = app.getTarget();
            H = ifftshift(r.otf);
            Y = real(ifft2(fft2(tg) .* H));
            Y = min(max(Y, 0), 1);
            if o.noise > 0
                Y = min(max(Y + o.noise * randn(size(Y)), 0), 1);
            end
            r.simT = tg;
            r.simY = Y;
            app.res = r;
        end

        % ================= render =================
        function render(app)
            %RENDER Update all graphics objects in place (never recreate).
            o = app.opts;
            r = app.res;
            c = app.cache;
            if isempty(r) || ~isfield(r, 'PSF')
                return;
            end
            N = o.N;
            cc = floor(N / 2) + 1;
            idx = (1:N) - cc;

            % ---- axis vectors (unit toggles) ----
            xUm = idx * r.pitchUm;
            if strcmp(o.psfUnits, 'um')
                xP = xUm; lblP = 'microns';
            else
                xP = idx / o.q; lblP = 'lambda*F#';
            end
            if strcmp(o.otfUnits, 'cmm')
                xF = idx * r.freqPitch; lblF = 'cycles/mm';
                xfac = r.freqPitch;
            else
                xF = idx * (o.q / N); lblF = 'f / cutoff';
                xfac = o.q / N;
            end
            xG = idx * (2 * o.q / N);   % pupil-plane grid units

            % ---- (1) pupil amplitude / phase ----
            set(app.h.imgPupil, 'XData', [xG(1) xG(end)], ...
                'YData', [xG(1) xG(end)]);
            if strcmp(o.pupilView, 'amp')
                set(app.h.imgPupil, 'CData', c.A);
                caxis(app.ax.pupil, [0 1]);
                app.ax.pupil.Title.String = 'Pupil amplitude';
            else
                ph = NaN(N);
                ph(c.mask) = r.W(c.mask);
                set(app.h.imgPupil, 'CData', ph);
                m = max(abs(r.W(c.mask)));
                if m < 1e-6, m = 0.01; end
                caxis(app.ax.pupil, [-m m]);
                app.ax.pupil.Title.String = 'Pupil phase W (waves)';
            end

            % ---- (2) PSF image: linear / log / gamma ----
            switch o.psfDisp
                case 'log'
                    D = (log10(max(r.PSF, 1e-6)) + 6) / 6;
                    psfTtl = 'PSF (log, 6 decades)';
                case 'gamma'
                    D = r.PSF .^ o.gamma;
                    psfTtl = sprintf('PSF (gamma %.2f)', o.gamma);
                otherwise
                    D = r.PSF;
                    psfTtl = 'PSF (linear)';
            end
            set(app.h.imgPSF, 'XData', [xP(1) xP(end)], ...
                'YData', [xP(1) xP(end)], 'CData', D);
            caxis(app.ax.psf, [0 1]);
            dP = xP(2) - xP(1);
            app.ax.psf.XLim = [xP(1) - dP / 2, xP(end) + dP / 2];
            app.ax.psf.YLim = app.ax.psf.XLim;
            app.ax.psf.XLabel.String = lblP;
            app.ax.psf.Title.String = psfTtl;

            % ---- (3) PSF cross-sections through the peak ----
            [~, pk] = max(r.PSF(:));
            [pr, pc] = ind2sub([N N], pk);
            set(app.h.sx, 'XData', xP, 'YData', r.PSF(pr, :)');
            set(app.h.sy, 'XData', xP, 'YData', r.PSF(:, pc)');
            [~, ip] = max(r.idealPSF(:));
            [ir, ic] = ind2sub([N N], ip);
            set(app.h.ix, 'XData', xP, 'YData', r.idealPSF(ir, :)');
            set(app.h.iy, 'XData', xP, 'YData', r.idealPSF(:, ic)');
            if o.showIdeal, visI = 'on'; else, visI = 'off'; end
            app.h.ix.Visible = visI;
            app.h.iy.Visible = visI;
            [visR, xr, xs, ys] = app.refVectors('psf');
            app.h.rx.Visible = visR;
            app.h.ry.Visible = visR;
            if strcmp(visR, 'on')
                set(app.h.rx, 'XData', xr, 'YData', xs);
                set(app.h.ry, 'XData', xr, 'YData', ys);
            end
            app.ax.psf1d.XLabel.String = lblP;
            app.ax.psf1d.YLim = [0 1.05];

            app.render2(xF, xfac, lblF, visI, visR, N, cc);
        end

        function [vis, xr, sx, sy] = refVectors(app, kind)
            %REFVECTORS Reference (frozen) curves in the CURRENT display
            % units; physical scales stored at freeze time are used, so a
            % frozen curve stays physically correct across lambda changes.
            rf = app.ref;
            o = app.opts;
            if isempty(rf)
                vis = 'off'; xr = NaN; sx = NaN; sy = NaN;
                return;
            end
            vis = 'on';
            i0 = (1:rf.N) - floor(rf.N / 2) - 1;
            switch kind
                case 'psf'
                    if strcmp(o.psfUnits, 'um')
                        xr = i0 * rf.pitchUm;
                    else
                        xr = i0 / rf.q;
                    end
                    sx = rf.sx;
                    sy = rf.sy;
                case 'mtf'
                    if strcmp(o.otfUnits, 'cmm')
                        xr = i0 * rf.freqPitch;
                    else
                        xr = i0 * (rf.q / rf.N);
                    end
                    sx = rf.sag;
                    sy = rf.tan;
                case 'ee'
                    if strcmp(o.psfUnits, 'um')
                        xr = rf.eeR' * rf.pitchUm;
                    else
                        xr = rf.eeR' / rf.q;
                    end
                    sx = rf.ee';
                    sy = [];
            end
            sx = sx(:)';
            if ~isempty(sy), sy = sy(:)'; end
        end

        function render2(app, xF, xfac, lblF, visI, visR, N, cc)
            %RENDER2 Second half of the render: OTF views, EE, simulation.
            o = app.opts;
            r = app.res;
            c = app.cache;

            % ---- (4) 2D MTF/PTF with cutoff circle ----
            if strcmp(o.otfView, 'mtf')
                D2 = r.mtf;
                caxis(app.ax.mtf2d, [0 1]);
                t2 = '2D MTF (cutoff circle)';
            else
                D2 = r.ptf / (2 * pi);      % waves, masked region = 0
                caxis(app.ax.mtf2d, [-0.5 0.5]);
                t2 = '2D PTF (waves)';
            end
            set(app.h.imgMTF, 'XData', [xF(1) xF(end)], ...
                'YData', [xF(1) xF(end)], 'CData', D2);
            dF = xF(2) - xF(1);
            app.ax.mtf2d.XLim = [xF(1) - dF / 2, xF(end) + dF / 2];
            app.ax.mtf2d.YLim = app.ax.mtf2d.XLim;
            if strcmp(o.otfUnits, 'cmm'), rad = r.cutoff; else, rad = 1; end
            th = 0:(2 * pi / 180):(2 * pi);
            set(app.h.circ, 'XData', rad * cos(th), 'YData', rad * sin(th));
            app.ax.mtf2d.XLabel.String = lblF;
            app.ax.mtf2d.Title.String = t2;

            % ---- (5) MTF cross-sections ----
            set(app.h.mSag, 'XData', xF, 'YData', r.mtf(cc, :)');
            set(app.h.mTan, 'XData', xF, 'YData', r.mtf(:, cc)');
            K = floor(N / 2) - 1;
            t = -K:K;
            sq = sqrt(0.5);
            set(app.h.mDia, 'XData', t * xfac, 'YData', ...
                interp2(r.mtf, cc + t * sq, cc + t * sq, 'linear'));
            iM = c.idealMtf;
            set(app.h.iSag, 'XData', xF, 'YData', iM(cc, :)');
            set(app.h.iTan, 'XData', xF, 'YData', iM(:, cc)');
            set(app.h.iDia, 'XData', t * xfac, 'YData', ...
                interp2(iM, cc + t * sq, cc + t * sq, 'linear'));
            app.h.iSag.Visible = visI;
            app.h.iTan.Visible = visI;
            app.h.iDia.Visible = visI;
            app.h.rSag.Visible = visR;
            app.h.rTan.Visible = visR;
            app.h.rDia.Visible = visR;
            if strcmp(visR, 'on')
                [~, xr, xs, ys] = app.refVectors('mtf');
                set(app.h.rSag, 'XData', xr, 'YData', xs);
                set(app.h.rTan, 'XData', xr, 'YData', ys);
                rf = app.ref;
                rK = floor(rf.N / 2) - 1;
                rt = -rK:rK;
                if strcmp(o.otfUnits, 'cmm')
                    rxf = rf.freqPitch;
                else
                    rxf = rf.q / rf.N;
                end
                set(app.h.rDia, 'XData', rt * rxf, 'YData', rf.diag(:)');
            end
            app.ax.mtf1d.XLabel.String = lblF;
            app.ax.mtf1d.YLim = [0 1.02];

            % ---- (6) encircled energy ----
            if strcmp(o.psfUnits, 'um')
                xE = r.eeR(:)' * r.pitchUm;
                lblE = 'radius (microns)';
            else
                xE = r.eeR(:)' / o.q;
                lblE = 'radius (lambda*F#)';
            end
            set(app.h.ee, 'XData', xE, 'YData', r.ee(:)');
            set(app.h.eeI, 'XData', xE, 'YData', c.idealEE(:)');
            app.h.eeI.Visible = visI;
            set(app.h.ee08, 'XData', [xE(1) xE(end)], 'YData', [0.8 0.8]);
            r80 = r.m.ee80Um / 2;                 % EE80 radius, microns
            if strcmp(o.psfUnits, 'um'), x80 = r80; else, x80 = r80 / r.scaleUm; end
            if isnan(x80)
                set(app.h.ee80, 'XData', NaN, 'YData', NaN);
            else
                set(app.h.ee80, 'XData', [x80 x80], 'YData', [0 1]);
            end
            app.h.eeR.Visible = visR;
            if strcmp(visR, 'on')
                [~, xrr, eer] = app.refVectors('ee');
                set(app.h.eeR, 'XData', xrr, 'YData', eer);
            end
            app.ax.ee.XLabel.String = lblE;
            app.ax.ee.YLim = [0 1.02];

            % ---- image simulation view ----
            set(app.h.imgTgt, 'XData', [1 N], 'YData', [1 N], ...
                'CData', r.simT);
            set(app.h.imgSim, 'XData', [1 N], 'YData', [1 N], ...
                'CData', r.simY);
            app.ax.tgt.XLim = [0.5 N + 0.5]; app.ax.tgt.YLim = [0.5 N + 0.5];
            app.ax.sim.XLim = [0.5 N + 0.5]; app.ax.sim.YLim = [0.5 N + 0.5];

            % ---- status, warnings, legends ----
            app.updateStatus();
            app.updateWarnings();
            app.updateLegends();
            if isempty(app.ref)
                app.ui.clearBtn.Enable = 'off';
            else
                app.ui.clearBtn.Enable = 'on';
            end
        end

        % ================= status / warnings / legends =================
        function updateStatus(app)
            m = app.res.m;
            app.StatusLbl.Text = sprintf( ...
                ['Strehl %.3f (Marechal %.3f)   |   RMS WFE %s waves ' ...
                 '(%s nm)   |   FWHM %s um (%s lF)   |   ' ...
                 'EE80 %s um (%s lF)   |   MTF50 %s c/mm (%s fc)   |   ' ...
                 'MTF10 %s c/mm (%s fc)'], ...
                m.strehl, m.marechal, ...
                app.fn(m.rmsWaves, '%.4f'), app.fn(m.rmsNm, '%.1f'), ...
                app.fn(m.fwhmUm, '%.3f'), app.fn(m.fwhmLf, '%.2f'), ...
                app.fn(m.ee80Um, '%.3f'), app.fn(m.ee80Lf, '%.2f'), ...
                app.fn(m.mtf50, '%.1f'), app.fn(m.mtf50Norm, '%.3f'), ...
                app.fn(m.mtf10, '%.1f'), app.fn(m.mtf10Norm, '%.3f'));
        end

        function s = fn(~, v, f)
            %FN Format a scalar, 'n/a' for NaN.
            if isnan(v), s = 'n/a'; else, s = sprintf(f, v); end
        end

        function updateWarnings(app)
            %UPDATEWARNINGS Non-blocking orange warnings (aliasing and
            % undersampling), shown in the status area - not dialogs.
            m = app.res.m;
            o = app.opts;
            w = {};
            if m.edgeFrac > 1e-3
                w{end+1} = sprintf(['PSF energy at grid edge (%.1e) - ' ...
                    'increase q to avoid wraparound'], m.edgeFrac);
            end
            if m.nyqRatio >= 0.9
                w{end+1} = sprintf(['MTF reaches the Nyquist limit of ' ...
                    'PSF sampling (q = %d, cutoff/Nyquist = %.2f)'], ...
                    o.q, m.nyqRatio);
            end
            if isempty(w)
                app.WarnLbl.Visible = 'off';
            else
                app.WarnLbl.Text = strjoin(w, '    |    ');
                app.WarnLbl.Visible = 'on';
            end
        end

        function updateLegends(app)
            %UPDATELEGENDS Rebuild the three 1-D legends when the overlay
            % set changes (ideal on/off, reference frozen/cleared).
            o = app.opts;
            hasRef = ~isempty(app.ref);
            key = [double(o.showIdeal), double(hasRef)];
            if isequal(app.legKey, key)
                return;
            end
            app.legKey = key;
            h = app.h;

            L = [h.sx h.sy]; N = {'x', 'y'};
            if o.showIdeal
                L = [L h.ix h.iy];
                N = [N {'ideal x', 'ideal y'}];
            end
            if hasRef
                L = [L h.rx h.ry];
                N = [N {'ref x', 'ref y'}];
            end
            legend(app.ax.psf1d, 'off');
            legend(app.ax.psf1d, L, N, 'Location', 'northeast', ...
                'FontSize', 8);

            L = [h.mSag h.mTan h.mDia];
            N = {'sagittal', 'tangential', 'diagonal'};
            if o.showIdeal
                L = [L h.iSag h.iTan h.iDia];
                N = [N {'ideal sag', 'ideal tan', 'ideal diag'}];
            end
            if hasRef
                L = [L h.rSag h.rTan h.rDia];
                N = [N {'ref sag', 'ref tan', 'ref diag'}];
            end
            legend(app.ax.mtf1d, 'off');
            legend(app.ax.mtf1d, L, N, 'Location', 'northeast', ...
                'FontSize', 8);

            L = h.ee; N = {'EE'};
            if o.showIdeal
                L = [L h.eeI]; N = [N {'ideal'}];
            end
            if hasRef
                L = [L h.eeR]; N = [N {'ref'}];
            end
            legend(app.ax.ee, 'off');
            legend(app.ax.ee, L, N, 'Location', 'southeast', 'FontSize', 8);
        end

        % ================= compare mode =================
        function onFreeze(app)
            %ONFREEZE Freeze the current 1-D curves as a dashed reference
            % overlay ("Freeze reference" button).
            r = app.res;
            if isempty(r) || ~isfield(r, 'PSF')
                return;
            end
            o = app.opts;
            N = o.N;
            cc = floor(N / 2) + 1;
            rf = struct();
            rf.N = N;
            rf.q = o.q;
            rf.pitchUm = r.pitchUm;
            rf.freqPitch = r.freqPitch;
            rf.cutoff = r.cutoff;
            [~, pk] = max(r.PSF(:));
            [pr, pc] = ind2sub([N N], pk);
            rf.sx = r.PSF(pr, :);
            rf.sy = r.PSF(:, pc);
            rf.sag = r.mtf(cc, :);
            rf.tan = r.mtf(:, cc);
            K = floor(N / 2) - 1;
            t = -K:K;
            sq = sqrt(0.5);
            rf.diag = interp2(r.mtf, cc + t * sq, cc + t * sq, 'linear');
            rf.ee = r.ee;
            rf.eeR = r.eeR;
            rf.label = app.refLabel();
            app.ref = rf;
            app.legKey = [];    % force legend rebuild on next render
            app.ax.psf1d.Title.String = sprintf( ...
                'PSF cross-sections  [ref: %s]', rf.label);
            app.markDirty();
        end

        function onClearRef(app)
            app.ref = [];
            app.legKey = [];
            app.ax.psf1d.Title.String = 'PSF cross-sections';
            app.markDirty();
        end

        function s = refLabel(app)
            %REFLABEL Short description of the frozen parameter set.
            o = app.opts;
            p = {o.shape};
            if strcmp(o.shape, 'annular')
                p{end+1} = sprintf('eps=%.2f', o.obstruction);
            end
            if ~strcmp(o.apod, 'uniform'), p{end+1} = o.apod; end
            if o.vanes > 0, p{end+1} = sprintf('%d vanes', o.vanes); end
            p{end+1} = sprintf('%.0f nm', o.lambda);
            p{end+1} = sprintf('f/%.1f', o.fnum);
            p{end+1} = sprintf('N=%d, q=%d', o.N, o.q);
            nz = find(o.coeff ~= 0);
            for k = 1:min(3, numel(nz))
                p{end+1} = sprintf('j%d=%+.3f', nz(k) - 1, o.coeff(nz(k)));
            end
            s = strjoin(p, ', ');
        end

        % ================= image-simulation target =================
        function img = getTarget(app)
            %GETTARGET Cached synthetic test target for the current N.
            tc = app.targetCache;
            if strcmp(tc.kind, app.opts.target) && tc.N == app.opts.N ...
                    && ~isempty(tc.img)
                img = tc.img;
                return;
            end
            img = psfx.makeTestTarget(app.opts.target, app.opts.N);
            app.targetCache = struct('kind', app.opts.target, ...
                'N', app.opts.N, 'img', img);
        end

        % ================= user mask =================
        function onLoadMask(app)
            %ONLOADMASK Load a user aperture mask from a .mat file or an
            % image file (PNG/JPG/TIF/BMP). Real I/O failures use uialert.
            [f, p] = uigetfile( ...
                {'*.mat;*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', ...
                 'Mask files (*.mat, images)'; '*.*', 'All files'}, ...
                'Load pupil mask');
            if isequal(f, 0)
                return;
            end
            fp = fullfile(p, f);
            try
                [~, ~, ext] = fileparts(f);
                if strcmpi(ext, '.mat')
                    S = load(fp);
                    fn = fieldnames(S);
                    M = [];
                    for k = 1:numel(fn)
                        v = S.(fn{k});
                        if (isnumeric(v) || islogical(v)) && ismatrix(v) ...
                                && ~isvector(v) && ~isempty(v)
                            M = double(v);
                            break;
                        end
                    end
                    if isempty(M)
                        error('PupilBench:mask', ...
                            'No 2-D numeric matrix found in the .mat file.');
                    end
                else
                    [img, map] = imread(fp);
                    if ~isempty(map)
                        img = ind2rgb(img, map);
                    end
                    if ndims(img) == 3
                        img = mean(double(img), 3);
                    else
                        img = double(img);
                    end
                    M = img;
                end
                m = max(M(:));
                if m > 1
                    M = M / m;          % accept 0..255 / 0..65535 images
                end
                M = min(max(M, 0), 1);
                app.userMaskRaw = M;
                app.userMaskN = [];
                app.userMaskId = app.userMaskId + 1;
                app.userMaskName = f;
                app.ui.maskLbl.Text = f;
                app.opts.shape = 'user';
                app.ui.shape.drop.Value = 'user';
                app.updateDeps();
                app.showMsg('');
                app.markDirty();
            catch ME
                uialert(app.Fig, ME.message, 'Could not load mask', ...
                    'Icon', 'error');
            end
        end

        % ================= exports =================
        function onExportPNG(app)
            [f, p] = uiputfile('*.png', 'Export figure PNG', ...
                'pupilbench.png');
            if isequal(f, 0), return; end
            fp = fullfile(p, f);
            ok = false;
            try
                exportgraphics(app.RightG, fp);   % container export
                ok = true;
            catch
            end
            if ~ok
                try
                    fr = getframe(app.Fig);       % fallback: rasterize
                    imwrite(fr.cdata, fp);
                    ok = true;
                catch
                end
            end
            if ~ok
                uialert(app.Fig, ['PNG export failed in this MATLAB ' ...
                    'version; try Export MAT / CSV instead.'], ...
                    'Export failed', 'Icon', 'error');
            end
        end

        function onExportMAT(app)
            [f, p] = uiputfile('*.mat', 'Export results struct', ...
                'pupilbench_results.mat');
            if isequal(f, 0), return; end
            try
                S = app.exportStruct();
                save(fullfile(p, f), '-struct', 'S');
            catch ME
                uialert(app.Fig, ME.message, 'Export failed', ...
                    'Icon', 'error');
            end
        end

        function onExportCSV(app)
            %ONEXPORTCSV MTF curves (sag/tan/diag + ideal) vs frequency.
            [f, p] = uiputfile('*.csv', 'Export MTF curves', ...
                'pupilbench_mtf.csv');
            if isequal(f, 0), return; end
            o = app.opts;
            r = app.res;
            c = app.cache;
            cc = floor(o.N / 2) + 1;
            idx = (1:o.N) - cc;
            K = floor(o.N / 2) - 1;
            t = -K:K;
            sq = sqrt(0.5);
            dg = interp2(r.mtf, cc + t * sq, cc + t * sq, 'linear');
            iM = c.idealMtf;
            idg = interp2(iM, cc + t * sq, cc + t * sq, 'linear');
            sel = abs(idx) <= K;
            dcol = NaN(o.N, 1);
            dcol(sel) = dg(idx(sel) + K + 1);
            icol = NaN(o.N, 1);
            icol(sel) = idg(idx(sel) + K + 1);
            sag = r.mtf(cc, :);
            tv = r.mtf(:, cc)';
            isag = iM(cc, :);
            itv = iM(:, cc)';
            fid = fopen(fullfile(p, f), 'w');
            if fid < 0
                uialert(app.Fig, 'Cannot open the output file.', ...
                    'Export failed', 'Icon', 'error');
                return;
            end
            try
                fprintf(fid, ['freq_norm,freq_cyc_per_mm,MTF_sag,' ...
                    'MTF_tan,MTF_diag,MTF_ideal_sag,MTF_ideal_tan,' ...
                    'MTF_ideal_diag\n']);
                for k = 1:o.N
                    fprintf(fid, ...
                        ['%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n'], ...
                        idx(k) * (o.q / o.N), idx(k) * r.freqPitch, ...
                        sag(k), tv(k), dcol(k), isag(k), itv(k), icol(k));
                end
            catch ME
                fclose(fid);
                uialert(app.Fig, ME.message, 'Export failed', ...
                    'Icon', 'error');
                return;
            end
            fclose(fid);
        end

















    end
end
