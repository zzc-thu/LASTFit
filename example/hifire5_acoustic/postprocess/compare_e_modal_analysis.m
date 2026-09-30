function compare_e_modal_analysis(mode)
% Fourier/POD/DMD comparison for E-2/E-3/E-4 wall pressure perturbations.

if nargin < 1
    mode = 'full';
end
spectraOnly = strcmpi(mode, 'spectra-only');

scriptDir = fileparts(mfilename('fullpath'));
caseReleaseDir = fileparts(scriptDir);
rootDir = getenv('LASTFIT_HIFIRE_RAW_DIR');
if isempty(rootDir)
    error(['Set LASTFIT_HIFIRE_RAW_DIR to the directory containing ', ...
        'E-2, E-3, E-4, and RESU (for example, Un_Hifire_fine).']);
end
outDir = fullfile(caseReleaseDir, 'output', 'manuscript_fine', ...
    'modal_comparison');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

cases = struct( ...
    'name', {'E-2', 'E-3', 'E-4'}, ...
    'epsilon', {5e-2, 5e-3, 5e-4}, ...
    'epsilonSource', {'inferred from folder name', 'inferred from folder name', 'Config.cfg'});

dtSaveDefault = 2e-4;
f0 = 666.6666666667;
targetKs = [41, 121, 1, 81];
axisNames = {'Long axis +', 'Long axis -', 'Short axis +', 'Short axis -'};
stationList = [40, 81, 121];

results = struct([]);
for ic = 1:numel(cases)
    fprintf('\n===== Modal analysis for %s =====\n', cases(ic).name);
    caseDir = fullfile(rootDir, cases(ic).name);
    cfgPath = fullfile(caseDir, 'Config.cfg');
    dtSave = dtSaveDefault;
    if exist(cfgPath, 'file')
        cfg = read_config(cfgPath);
        cases(ic).epsilon = get_cfg_number(cfg, 'epsilon', cases(ic).epsilon);
        cases(ic).epsilonSource = 'Config.cfg';
        dt = get_cfg_number(cfg, 'FixedDeltaTime', 1e-7);
        steps = get_cfg_number(cfg, 'StepsWriteData', 2000);
        dtSave = dt * steps;
    end

    result = load_case_meridians(caseDir, cases(ic), dtSave, targetKs, axisNames);
    result = analyze_fourier(result, f0, stationList);
    if ~spectraOnly
        result = analyze_pod(result);
        result = analyze_dmd(result, f0);
    end
    results = [results; result]; %#ok<AGROW>
end

if spectraOnly
    plot_station_spectra(outDir, results, f0, stationList);
    compareDir = fullfile(rootDir, 'Compare');
    if ~exist(compareDir, 'dir')
        mkdir(compareDir);
    end
    copyfile(fullfile(outDir, 'fourier_station_spectra_by_axis.png'), ...
        fullfile(compareDir, 'fourier_station_spectra_by_axis.png'), 'f');
    fprintf('\nDone. Updated spectra figure only:\n%s\n', ...
        fullfile(compareDir, 'fourier_station_spectra_by_axis.png'));
    return;
end

write_metrics_csv(outDir, results);
write_station_metrics_csv(outDir, results);
plot_epsilon_summary(outDir, results);
plot_station_spectra(outDir, results, f0, stationList);
plot_station_amplitude_bars(outDir, results, stationList);
plot_station_harmonic_ratio_bars(outDir, results, stationList);
plot_f0_profiles(outDir, results);
plot_h2_profiles(outDir, results);
plot_pod_energy(outDir, results);
plot_pod_modes(outDir, results);
plot_pod_time_coefficients(outDir, results);
plot_dmd_eigs(outDir, results);
plot_dmd_frequency_growth(outDir, results);
plot_dmd_reconstruction_error(outDir, results);
plot_dmd_f0_modes(outDir, results);
write_report_v2(outDir, results, f0, stationList);
sync_compare_outputs(rootDir, outDir, results, f0, stationList);

fprintf('\nDone. Output directory:\n%s\n', outDir);
end

function result = load_case_meridians(caseDir, caseInfo, dtSave, targetKs, axisNames)
pertDir = fullfile(caseDir, 'Pert');
masters = dir(fullfile(pertDir, 'Pert_Result_*.pvts'));
if isempty(masters)
    masters = dir(fullfile(pertDir, 'Result_*.pvts'));
end
if isempty(masters)
    error('No Pert PVTS files found in %s.', pertDir);
end
[frameIds, order] = sort(extract_frame_ids({masters.name}));
masters = masters(order);
nFrames = numel(masters);
t = (0:nFrames-1)' * dtSave;
lateIdx = ceil(nFrames/2):nFrames;

nAxes = numel(targetKs);
pAxes = cell(nAxes, 1);
iAxes = cell(nAxes, 1);
xAxes = cell(nAxes, 1);

for it = 1:nFrames
    masterPath = fullfile(pertDir, masters(it).name);
    needPoints = it == 1;
    axesFrame = extract_wall_axes(masterPath, targetKs, needPoints);
    for ia = 1:nAxes
        if it == 1
            iAxes{ia} = axesFrame(ia).i(:);
            xAxes{ia} = axesFrame(ia).x(:);
            pAxes{ia} = nan(numel(iAxes{ia}), nFrames);
        end
        [tf, loc] = ismember(iAxes{ia}, axesFrame(ia).i(:));
        p = nan(size(iAxes{ia}));
        p(tf) = axesFrame(ia).p(loc(tf));
        pAxes{ia}(:, it) = p;
    end
    if mod(it, max(1, floor(nFrames/8))) == 0 || it == nFrames
        fprintf('%s frame %d / %d\n', caseInfo.name, it, nFrames);
    end
end

[Xfull, meta] = stack_axes(pAxes, iAxes, xAxes, axisNames);
Xfull = fill_nan_matrix(Xfull);
Xlate = Xfull(:, lateIdx);
Xlate = Xlate - mean(Xlate, 2);

result.name = caseInfo.name;
result.epsilon = caseInfo.epsilon;
result.epsilonSource = caseInfo.epsilonSource;
result.caseDir = caseDir;
result.frameIds = frameIds;
result.nFrames = nFrames;
result.dtSave = dtSave;
result.fs = 1 / dtSave;
result.t = t;
result.lateIdx = lateIdx(:);
result.tLate = t(lateIdx);
result.targetKs = targetKs;
result.axisNames = axisNames;
result.pAxes = pAxes;
result.iAxes = iAxes;
result.xAxes = xAxes;
result.Xfull = Xfull;
result.Xlate = Xlate;
result.meta = meta;
end

function result = analyze_fourier(result, f0, stationList)
f2 = 2 * f0;
tLate = result.tLate(:);
nLate = numel(tLate);
phase1 = exp(-1i * 2*pi*f0*tLate);
phase2 = exp(-1i * 2*pi*f2*tLate);

nAxes = numel(result.pAxes);
amp1 = cell(nAxes, 1);
amp2 = cell(nAxes, 1);
phase1Rad = cell(nAxes, 1);
rmsLate = cell(nAxes, 1);
for ia = 1:nAxes
    y = result.pAxes{ia}(:, result.lateIdx);
    y = fill_nan_matrix(y);
    y = y - mean(y, 2);
    coeff1 = (2 / nLate) * (y * phase1);
    coeff2 = (2 / nLate) * (y * phase2);
    amp1{ia} = abs(coeff1);
    amp2{ia} = abs(coeff2);
    phase1Rad{ia} = unwrap(angle(coeff1));
    rmsLate{ia} = sqrt(mean(y.^2, 2));
end

longIdx = 1:2;
shortIdx = 3:4;
longAmp = mean_same_grid(amp1(longIdx));
shortAmp = mean_same_grid(amp1(shortIdx));
longAmp2 = mean_same_grid(amp2(longIdx));
shortAmp2 = mean_same_grid(amp2(shortIdx));
longRms = mean_same_grid(rmsLate(longIdx));
shortRms = mean_same_grid(rmsLate(shortIdx));
longPhase = circular_mean(phase1Rad(longIdx));
shortPhase = circular_mean(phase1Rad(shortIdx));

iCommon = result.iAxes{1};
xCommon = result.xAxes{1};
validMask = iCommon >= 10;

station = station_fourier_metrics(result.pAxes, result.iAxes, result.lateIdx, result.t, result.fs, stationList, f0);
midStation = stationList(ceil(numel(stationList) / 2));
[freqLong, specLong, rawLong] = station_spectrum(result.pAxes, result.iAxes, result.lateIdx, result.t, result.fs, midStation, 1:2);
[freqShort, specShort, rawShort] = station_spectrum(result.pAxes, result.iAxes, result.lateIdx, result.t, result.fs, midStation, 3:4);

result.fourier.amp1 = amp1;
result.fourier.amp2 = amp2;
result.fourier.phase1Rad = phase1Rad;
result.fourier.rmsLate = rmsLate;
result.fourier.i = iCommon;
result.fourier.x = xCommon;
result.fourier.longAmp = longAmp;
result.fourier.shortAmp = shortAmp;
result.fourier.longAmp2 = longAmp2;
result.fourier.shortAmp2 = shortAmp2;
result.fourier.longRms = longRms;
result.fourier.shortRms = shortRms;
result.fourier.longPhase = longPhase;
result.fourier.shortPhase = shortPhase;
result.fourier.meanLongAmp = mean(longAmp(validMask), 'omitnan');
result.fourier.meanShortAmp = mean(shortAmp(validMask), 'omitnan');
result.fourier.meanLongAmp2 = mean(longAmp2(validMask), 'omitnan');
result.fourier.meanShortAmp2 = mean(shortAmp2(validMask), 'omitnan');
result.fourier.meanLongRms = mean(longRms(validMask), 'omitnan');
result.fourier.meanShortRms = mean(shortRms(validMask), 'omitnan');
result.fourier.longShortAmpRatio = result.fourier.meanLongAmp / result.fourier.meanShortAmp;
result.fourier.longH2H1 = result.fourier.meanLongAmp2 / result.fourier.meanLongAmp;
result.fourier.shortH2H1 = result.fourier.meanShortAmp2 / result.fourier.meanShortAmp;
result.fourier.freqLong = freqLong;
result.fourier.specLong = specLong;
result.fourier.rawLong = rawLong;
result.fourier.freqShort = freqShort;
result.fourier.specShort = specShort;
result.fourier.rawShort = rawShort;
result.fourier.station = station;
end

function result = analyze_pod(result)
X = result.Xlate;
[U, S, V] = svd(X, 'econ');
s = diag(S);
energy = s.^2 / sum(s.^2);
result.pod.U = U;
result.pod.S = S;
result.pod.V = V;
result.pod.singularValues = s;
result.pod.energy = energy;
result.pod.cumulativeEnergy = cumsum(energy);
result.pod.coeff = S * V';
result.pod.n90 = find(result.pod.cumulativeEnergy >= 0.90, 1, 'first');
result.pod.n95 = find(result.pod.cumulativeEnergy >= 0.95, 1, 'first');
end

function result = analyze_dmd(result, f0)
X = result.Xlate;
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);
[U, S, V] = svd(X1, 'econ');
s = diag(S);
tol = max(size(X1)) * eps(max(s));
rAuto = sum(s > tol);
r = min([12, rAuto, size(X1, 2)]);
if r < 2
    r = min(2, size(X1, 2));
end
U = U(:, 1:r);
S = S(1:r, 1:r);
V = V(:, 1:r);
Atilde = U' * X2 * V / S;
[W, D] = eig(Atilde);
lambda = diag(D);
Phi = X2 * V / S * W;
omega = log(lambda) / result.dtSave;
freq = imag(omega) / (2*pi);
growth = real(omega);
b = Phi \ X(:, 1);

timeSteps = 0:(size(X, 2)-1);
dynamics = zeros(numel(lambda), numel(timeSteps));
for k = 1:numel(timeSteps)
    dynamics(:, k) = b .* (lambda .^ timeSteps(k));
end
Xrec = real(Phi * dynamics);
relError = norm(X - Xrec, 'fro') / max(norm(X, 'fro'), eps);

positive = find(freq >= 0);
[~, loc] = min(abs(freq(positive) - f0));
f0Mode = positive(loc);
[~, ampOrder] = sort(abs(b), 'descend');

result.dmd.rank = r;
result.dmd.lambda = lambda;
result.dmd.omega = omega;
result.dmd.freq = freq;
result.dmd.growth = growth;
result.dmd.Phi = Phi;
result.dmd.b = b;
result.dmd.relError = relError;
result.dmd.f0Mode = f0Mode;
result.dmd.f0ModeFrequency = freq(f0Mode);
result.dmd.f0ModeGrowth = growth(f0Mode);
result.dmd.amplitudeOrder = ampOrder;
end

function write_metrics_csv(outDir, results)
fn = fullfile(outDir, 'modal_comparison_metrics.csv');
fid = fopen(fn, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['case,epsilon,epsilon_source,n_frames,late_frame_count,dt_save,fs,' ...
    'mean_long_amp_f0,mean_short_amp_f0,long_short_amp_ratio,long_amp_over_epsilon,short_amp_over_epsilon,' ...
    'mean_long_amp_2f0,mean_short_amp_2f0,long_2f0_over_f0,short_2f0_over_f0,' ...
    'pod_energy_mode1,pod_energy_mode2,pod_modes_90,pod_modes_95,' ...
    'dmd_rank,dmd_f0_mode_frequency,dmd_f0_mode_growth,dmd_reconstruction_error\n']);
for ic = 1:numel(results)
    r = results(ic);
    fprintf(fid, '%s,%.12g,%s,%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d,%d,%d,%.12g,%.12g,%.12g\n', ...
        r.name, r.epsilon, r.epsilonSource, r.nFrames, numel(r.lateIdx), r.dtSave, r.fs, ...
        r.fourier.meanLongAmp, r.fourier.meanShortAmp, r.fourier.longShortAmpRatio, ...
        r.fourier.meanLongAmp / r.epsilon, r.fourier.meanShortAmp / r.epsilon, ...
        r.fourier.meanLongAmp2, r.fourier.meanShortAmp2, r.fourier.longH2H1, r.fourier.shortH2H1, ...
        r.pod.energy(1), get_or_nan(r.pod.energy, 2), r.pod.n90, r.pod.n95, ...
        r.dmd.rank, r.dmd.f0ModeFrequency, r.dmd.f0ModeGrowth, r.dmd.relError);
end
end

function write_station_metrics_csv(outDir, results)
fn = fullfile(outDir, 'station_fourier_metrics.csv');
fid = fopen(fn, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'case,epsilon,station_i,axis,amp_f0,amp_f0_over_epsilon,amp_2f0,amp_2f0_over_epsilon,h2_over_h1\n');
for ic = 1:numel(results)
    r = results(ic);
    for is = 1:numel(r.fourier.station)
        s = r.fourier.station(is);
        fprintf(fid, '%s,%.12g,%d,long,%.12g,%.12g,%.12g,%.12g,%.12g\n', ...
            r.name, r.epsilon, s.i, s.longAmpF0, s.longAmpF0 / r.epsilon, ...
            s.longAmp2F0, s.longAmp2F0 / r.epsilon, s.longH2H1);
        fprintf(fid, '%s,%.12g,%d,short,%.12g,%.12g,%.12g,%.12g,%.12g\n', ...
            r.name, r.epsilon, s.i, s.shortAmpF0, s.shortAmpF0 / r.epsilon, ...
            s.shortAmp2F0, s.shortAmp2F0 / r.epsilon, s.shortH2H1);
    end
end
end

function plot_epsilon_summary(outDir, results)
fig = make_figure([1500 900]);
epsVals = [results.epsilon];
longAmp = arrayfun(@(r) r.fourier.meanLongAmp, results);
shortAmp = arrayfun(@(r) r.fourier.meanShortAmp, results);
longH = arrayfun(@(r) r.fourier.longH2H1, results);
shortH = arrayfun(@(r) r.fourier.shortH2H1, results);
labels = case_labels(results);

subplot(1, 3, 1);
loglog(epsVals, longAmp, '-o', 'LineWidth', 2.2, 'MarkerSize', 8); hold on;
loglog(epsVals, shortAmp, '-s', 'LineWidth', 2.2, 'MarkerSize', 8);
grid on; xlabel('\epsilon'); ylabel('mean |P''_{f0}|');
title('Fundamental amplitude scaling');
legend({'Long axis true amplitude', 'Short axis true amplitude'}, 'Location', 'best');
style_pub(gca);
annotate_points(epsVals, longAmp, longAmp);
annotate_points(epsVals, shortAmp, shortAmp);

subplot(1, 3, 2);
semilogx(epsVals, longAmp ./ epsVals, '-o', 'LineWidth', 2.2, 'MarkerSize', 8); hold on;
semilogx(epsVals, shortAmp ./ epsVals, '-s', 'LineWidth', 2.2, 'MarkerSize', 8);
grid on; xlabel('\epsilon'); ylabel('mean |P''_{f0}| / \epsilon');
title('Normalized fundamental response');
legend({'Long axis', 'Short axis'}, 'Location', 'best');
style_pub(gca);
for i = 1:numel(results)
    text(epsVals(i), longAmp(i) ./ epsVals(i), labels{i}, 'FontSize', 10, 'VerticalAlignment', 'bottom');
end

subplot(1, 3, 3);
semilogx(epsVals, longH, '-o', 'LineWidth', 2.2, 'MarkerSize', 8); hold on;
semilogx(epsVals, shortH, '-s', 'LineWidth', 2.2, 'MarkerSize', 8);
grid on; xlabel('\epsilon'); ylabel('|P''_{2f0}| / |P''_{f0}|');
title('Second harmonic ratio');
legend({'Long axis', 'Short axis'}, 'Location', 'best');
style_pub(gca);
save_fig(fig, fullfile(outDir, 'epsilon_fourier_summary.png'));
end

function plot_station_spectra(outDir, results, f0, stationList)
fig = make_figure([2400 1700]);
colors = lines(numel(results));
for is = 1:numel(stationList)
    yMaxPair = station_spectrum_ymax(results, is);
    for axGroup = 1:2
        subplot(numel(stationList), 2, (is-1)*2 + axGroup);
        amps = zeros(numel(results), 1);
        for ic = 1:numel(results)
            s = results(ic).fourier.station(is);
            if axGroup == 1
                freq = s.freqLong;
                spec = s.specLong;
                ampF0 = s.longAmpF0;
                axisText = 'long axis';
            else
                freq = s.freqShort;
                spec = s.specShort;
                ampF0 = s.shortAmpF0;
                axisText = 'short axis';
            end
            plot(freq, spec ./ results(ic).epsilon, '-', 'LineWidth', 2.8, 'Color', colors(ic, :)); hold on;
            amps(ic) = ampF0;
        end
        xline(f0, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
        xline(2*f0, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, 'HandleVisibility', 'off');
        grid on; xlabel('frequency'); ylabel('amplitude / \epsilon');
        title(sprintf('i = %d, %s', stationList(is), axisText));
        xlim([0, 2300]);
        ylim([0, yMaxPair]);
        ampText = sprintf('|P''_{f0}| true: %.3g / %.3g / %.3g', amps(1), amps(2), amps(3));
        text(0.985, 0.90, ampText, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'normal', ...
            'BackgroundColor', 'w', 'Margin', 2);
        style_pub(gca);
    end
end
legendAx = axes('Position', [0.17 0.955 0.66 0.035], 'Visible', 'off');
hold(legendAx, 'on');
h = gobjects(numel(results), 1);
for ic = 1:numel(results)
    h(ic) = plot(legendAx, nan, nan, '-', 'LineWidth', 3.2, 'Color', colors(ic, :));
end
lgd = legend(legendAx, h, case_names(results), ...
    'Orientation', 'horizontal', 'Location', 'north', 'Box', 'off');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
save_fig(fig, fullfile(outDir, 'fourier_station_spectra_by_axis.png'));
end

function plot_station_amplitude_bars(outDir, results, stationList)
fig = make_figure([2100 1050]);
allVals = zeros(numel(stationList), numel(results), 2);
for ic = 1:numel(results)
    for is = 1:numel(stationList)
        s = results(ic).fourier.station(is);
        allVals(is, ic, 1) = s.longAmpF0;
        allVals(is, ic, 2) = s.shortAmpF0;
    end
end
yMaxCommon = max(allVals(:));
for axGroup = 1:2
    subplot(1, 2, axGroup);
    vals = zeros(numel(stationList), numel(results));
    for ic = 1:numel(results)
        for is = 1:numel(stationList)
            s = results(ic).fourier.station(is);
            if axGroup == 1
                vals(is, ic) = s.longAmpF0;
            else
                vals(is, ic) = s.shortAmpF0;
            end
        end
    end
    b = bar(vals, 'grouped'); %#ok<NASGU>
    set(gca, 'XTick', 1:numel(stationList), 'XTickLabel', compose('i=%d', stationList));
    grid on; xlabel('station'); ylabel('true |P''_{f0}|');
    if axGroup == 1
        title('Long-axis true |P''_{f0}|');
    else
        title('Short-axis true |P''_{f0}|');
    end
    lgd = legend(case_names(results), 'Location', 'northoutside', 'Orientation', 'horizontal');
    set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
    style_pub(gca);
    annotate_grouped_bars(vals, yMaxCommon);
end
save_fig(fig, fullfile(outDir, 'fourier_station_true_f0_amplitudes.png'));
end

function plot_station_harmonic_ratio_bars(outDir, results, stationList)
fig = make_figure([2100 1050]);
allVals = zeros(numel(stationList), numel(results), 2);
for ic = 1:numel(results)
    for is = 1:numel(stationList)
        s = results(ic).fourier.station(is);
        allVals(is, ic, 1) = s.longH2H1;
        allVals(is, ic, 2) = s.shortH2H1;
    end
end
yMaxCommon = max(allVals(:));
for axGroup = 1:2
    subplot(1, 2, axGroup);
    vals = zeros(numel(stationList), numel(results));
    amp2 = zeros(numel(stationList), numel(results));
    for ic = 1:numel(results)
        for is = 1:numel(stationList)
            s = results(ic).fourier.station(is);
            if axGroup == 1
                vals(is, ic) = s.longH2H1;
                amp2(is, ic) = s.longAmp2F0;
            else
                vals(is, ic) = s.shortH2H1;
                amp2(is, ic) = s.shortAmp2F0;
            end
        end
    end
    bar(vals, 'grouped');
    set(gca, 'XTick', 1:numel(stationList), 'XTickLabel', compose('i=%d', stationList));
    grid on; xlabel('station'); ylabel('|P''_{2f0}| / |P''_{f0}|');
    if axGroup == 1
        title('Long-axis |P''_{2f0}| / |P''_{f0}|');
    else
        title('Short-axis |P''_{2f0}| / |P''_{f0}|');
    end
    lgd = legend(case_names(results), 'Location', 'northoutside', 'Orientation', 'horizontal');
    set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
    style_pub(gca);
    annotate_grouped_bars(vals, yMaxCommon);
    text(0.02, 0.96, sprintf('Bar labels show H2/H1; true |P_2f0| range %.3g to %.3g', min(amp2(:)), max(amp2(:))), ...
        'Units', 'normalized', 'FontSize', 11, 'VerticalAlignment', 'top', 'BackgroundColor', 'w');
end
save_fig(fig, fullfile(outDir, 'fourier_station_second_harmonic_ratios.png'));
end

function plot_f0_profiles(outDir, results)
fig = make_figure([1500 1050]);
colors = lines(numel(results));
labels = case_labels(results);
yMaxTrue = axis_f0_ymax(results, false);
yMaxNorm = axis_f0_ymax(results, true);
subplot(2, 2, 1);
for ic = 1:numel(results)
    plot(results(ic).fourier.x, results(ic).fourier.longAmp, ...
        '-', 'LineWidth', 2.0, 'Color', colors(ic, :)); hold on;
end
grid on; xlabel('X'); ylabel('true |P''_{f0}|');
title('Long-axis fundamental true amplitude');
legend(labels, 'Location', 'best');
ylim([0, yMaxTrue]);
style_pub(gca);

subplot(2, 2, 2);
for ic = 1:numel(results)
    plot(results(ic).fourier.x, results(ic).fourier.longAmp ./ results(ic).epsilon, ...
        '-', 'LineWidth', 2.0, 'Color', colors(ic, :)); hold on;
end
grid on; xlabel('X'); ylabel('|P''_{f0}| / \epsilon');
title('Long-axis normalized fundamental amplitude');
legend(labels, 'Location', 'best');
ylim([0, yMaxNorm]);
style_pub(gca);

subplot(2, 2, 3);
for ic = 1:numel(results)
    plot(results(ic).fourier.x, results(ic).fourier.shortAmp, ...
        '-', 'LineWidth', 2.0, 'Color', colors(ic, :)); hold on;
end
grid on; xlabel('X'); ylabel('true |P''_{f0}|');
title('Short-axis fundamental true amplitude');
legend(labels, 'Location', 'best');
ylim([0, yMaxTrue]);
style_pub(gca);

subplot(2, 2, 4);
for ic = 1:numel(results)
    plot(results(ic).fourier.x, results(ic).fourier.shortAmp ./ results(ic).epsilon, ...
        '-', 'LineWidth', 2.0, 'Color', colors(ic, :)); hold on;
end
grid on; xlabel('X'); ylabel('|P''_{f0}| / \epsilon');
title('Short-axis normalized fundamental amplitude');
legend(labels, 'Location', 'best');
ylim([0, yMaxNorm]);
style_pub(gca);
save_fig(fig, fullfile(outDir, 'fourier_axis_f0_profiles.png'));
end

function plot_h2_profiles(outDir, results)
fig = make_figure([1500 1050]);
colors = lines(numel(results));
labels = case_labels(results);
yMaxTrue = axis_h2_ymax(results, false);
yMaxRatio = axis_h2_ymax(results, true);
subplot(2, 2, 1);
for ic = 1:numel(results)
    plot(results(ic).fourier.x, results(ic).fourier.longAmp2, '-', 'LineWidth', 2.0, 'Color', colors(ic, :)); hold on;
end
grid on; xlabel('X'); ylabel('true |P''_{2f0}|');
title('Long-axis second harmonic true amplitude');
legend(labels, 'Location', 'best');
ylim([0, yMaxTrue]);
style_pub(gca);

subplot(2, 2, 2);
for ic = 1:numel(results)
    ratio = results(ic).fourier.longAmp2 ./ max(results(ic).fourier.longAmp, eps);
    plot(results(ic).fourier.x, ratio, '-', 'LineWidth', 2.0, 'Color', colors(ic, :)); hold on;
end
grid on; xlabel('X'); ylabel('|P''_{2f0}| / |P''_{f0}|');
title('Long-axis second harmonic ratio');
legend(labels, 'Location', 'best');
ylim([0, yMaxRatio]);
style_pub(gca);

subplot(2, 2, 3);
for ic = 1:numel(results)
    plot(results(ic).fourier.x, results(ic).fourier.shortAmp2, '-', 'LineWidth', 2.0, 'Color', colors(ic, :)); hold on;
end
grid on; xlabel('X'); ylabel('true |P''_{2f0}|');
title('Short-axis second harmonic true amplitude');
legend(labels, 'Location', 'best');
ylim([0, yMaxTrue]);
style_pub(gca);

subplot(2, 2, 4);
for ic = 1:numel(results)
    ratio = results(ic).fourier.shortAmp2 ./ max(results(ic).fourier.shortAmp, eps);
    plot(results(ic).fourier.x, ratio, '-', 'LineWidth', 2.0, 'Color', colors(ic, :)); hold on;
end
grid on; xlabel('X'); ylabel('|P''_{2f0}| / |P''_{f0}|');
title('Short-axis second harmonic ratio');
legend(labels, 'Location', 'best');
ylim([0, yMaxRatio]);
style_pub(gca);
save_fig(fig, fullfile(outDir, 'fourier_axis_second_harmonic_ratio.png'));
end

function plot_pod_energy(outDir, results)
fig = make_figure([1250 780]);
colors = lines(numel(results));
labels = case_labels(results);
for ic = 1:numel(results)
    n = min(12, numel(results(ic).pod.energy));
    plot(1:n, results(ic).pod.cumulativeEnergy(1:n), '-o', ...
        'LineWidth', 2.0, 'MarkerSize', 6, 'Color', colors(ic, :)); hold on;
end
yline(0.90, 'k--', 'LineWidth', 1.2);
yline(0.95, 'k:', 'LineWidth', 1.2);
grid on; xlabel('POD mode count'); ylabel('cumulative energy');
title('POD cumulative energy in late window');
legend(labels, 'Location', 'southeast');
ylim([0, 1.02]);
style_pub(gca);
save_fig(fig, fullfile(outDir, 'pod_cumulative_energy.png'));
end

function plot_pod_modes(outDir, results)
fig = make_figure([2100 1250]);
colors = lines(numel(results));
labels = case_labels(results);
for m = 1:3
    subplot(3, 1, m);
    for ic = 1:numel(results)
        [x, y] = axis_component(results(ic).pod.U(:, m), results(ic).meta, 1);
        y = y / max(abs(y) + eps);
        plot(x, y, '-', 'LineWidth', 2.0, 'Color', colors(ic, :)); hold on;
    end
    grid on; xlabel('X'); ylabel(sprintf('mode %d', m));
    title(sprintf('POD mode %d shape on long axis +', m));
    if m == 1
        lgd = legend(labels, 'Location', 'northeastoutside');
        set(lgd, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
    end
    style_pub(gca);
end
save_fig(fig, fullfile(outDir, 'pod_long_axis_mode_shapes.png'));
end

function plot_pod_time_coefficients(outDir, results)
fig = make_figure([1500 900]);
colors = lines(numel(results));
labels = case_labels(results);
for m = 1:2
    subplot(2, 1, m);
    for ic = 1:numel(results)
        tt = results(ic).tLate - results(ic).tLate(1);
        a = results(ic).pod.coeff(m, :);
        a = a / max(abs(a) + eps);
        plot(tt, a, '-o', 'LineWidth', 1.8, 'MarkerSize', 4, 'Color', colors(ic, :)); hold on;
    end
    grid on; xlabel('t - t_{late,start}'); ylabel(sprintf('a_%d / max|a_%d|', m, m));
    title(sprintf('POD mode %d normalized time coefficient', m));
    legend(labels, 'Location', 'best');
    style_pub(gca);
end
save_fig(fig, fullfile(outDir, 'pod_time_coefficients.png'));
end

function plot_dmd_eigs(outDir, results)
fig = make_figure([900 850]);
theta = linspace(0, 2*pi, 400);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1.2); hold on;
colors = lines(numel(results));
labels = case_labels(results);
for ic = 1:numel(results)
    lam = results(ic).dmd.lambda;
    scatter(real(lam), imag(lam), 70, colors(ic, :), 'filled');
end
axis equal; grid on; xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('DMD eigenvalues');
legend([{'unit circle'}, labels], 'Location', 'bestoutside');
style_pub(gca);
save_fig(fig, fullfile(outDir, 'dmd_eigenvalues.png'));
end

function plot_dmd_frequency_growth(outDir, results)
fig = make_figure([1250 780]);
colors = lines(numel(results));
labels = case_labels(results);
for ic = 1:numel(results)
    freq = abs(results(ic).dmd.freq);
    growth = results(ic).dmd.growth;
    amp = abs(results(ic).dmd.b);
    amp = 35 + 130 * amp / max(amp + eps);
    scatter(freq, growth, amp, colors(ic, :), 'filled', 'MarkerFaceAlpha', 0.65); hold on;
end
grid on; xlabel('|frequency|'); ylabel('growth rate');
title('DMD frequency-growth map');
legend(labels, 'Location', 'best');
xlim([0, 2500]);
style_pub(gca);
save_fig(fig, fullfile(outDir, 'dmd_frequency_growth.png'));
end

function plot_dmd_reconstruction_error(outDir, results)
fig = make_figure([950 700]);
err = arrayfun(@(r) r.dmd.relError, results);
bar(err, 0.62);
set(gca, 'XTick', 1:numel(results), 'XTickLabel', case_labels(results));
grid on; ylabel('relative Frobenius error');
title('DMD reconstruction error in late window');
style_pub(gca);
save_fig(fig, fullfile(outDir, 'dmd_reconstruction_error.png'));
end

function plot_dmd_f0_modes(outDir, results)
fig = make_figure([2100 1250]);
colors = lines(numel(results));
labels = case_labels(results);
for ic = 1:numel(results)
    subplot(numel(results), 1, ic);
    modeVec = real(results(ic).dmd.Phi(:, results(ic).dmd.f0Mode));
    [x, y] = axis_component(modeVec, results(ic).meta, 1);
    y = y / max(abs(y) + eps);
    plot(x, y, '-', 'LineWidth', 2.2, 'Color', colors(ic, :));
    grid on; xlabel('X'); ylabel('mode');
    title(sprintf('%s DMD f0-nearest mode: %.3f Hz, growth %.3g', ...
        labels{ic}, results(ic).dmd.f0ModeFrequency, results(ic).dmd.f0ModeGrowth));
    style_pub(gca);
end
save_fig(fig, fullfile(outDir, 'dmd_f0_mode_shapes_long_axis.png'));
end

function write_report_v2(outDir, results, f0, stationList)
fn = fullfile(outDir, 'figure_analysis_report.md');
fid = fopen_utf8(fn);
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# 不同扰动幅值的 Fourier/POD/DMD 对比分析\n\n');
fprintf(fid, '本次后处理比较 `E-2`、`E-3`、`E-4` 三组快声波扰动幅值算例。分析对象为壁面 `j=1` 上四条代表母线：长轴正负向 `k=41,121` 与短轴正负向 `k=1,81` 的压力扰动 `P_pert`。为了避免初始瞬态主导模态，Fourier、POD 与 DMD 均使用每个算例后半段时间窗。\n\n');
fprintf(fid, '采样间隔为 `Delta t_save = %.6g`，目标基频为 `f0 = %.6g Hz`。Fourier 频谱和站位幅值分析覆盖 `i = %s` 三个站位，并分别比较长轴与短轴。\n\n', results(1).dtSave, f0, join_ints(stationList));
fprintf(fid, '为方便从图上直接比较长轴与短轴，所有包含长短轴并列或上下对比的图均使用成对一致的纵轴范围：频谱图每个站位的长/短轴共用 y 轴，站位柱状图的左右面板共用 y 轴，沿流向剖面图的长/短轴真实幅值面板和归一化面板分别共用 y 轴。\n\n');

fprintf(fid, '## 数据与主要指标\n\n');
fprintf(fid, '| case | epsilon | frames | late frames | long P_f0/epsilon | short P_f0/epsilon | long H2/H1 | short H2/H1 | POD E1 | POD modes 95%% | DMD error |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ic = 1:numel(results)
    r = results(ic);
    fprintf(fid, '| %s | %.3g | %d | %d | %.4g | %.4g | %.4g | %.4g | %.4g | %d | %.4g |\n', ...
        r.name, r.epsilon, r.nFrames, numel(r.lateIdx), ...
        r.fourier.meanLongAmp / r.epsilon, r.fourier.meanShortAmp / r.epsilon, ...
        r.fourier.longH2H1, r.fourier.shortH2H1, r.pod.energy(1), r.pod.n95, r.dmd.relError);
end

fprintf(fid, '\n下表给出各站位的真实主频幅值 `|P_f0|` 和二倍频比例。这里的 `|P_f0|` 未除以 `epsilon`，可直接反映实际压力扰动量级。\n\n');
fprintf(fid, '| case | station | long P_f0 | short P_f0 | long H2/H1 | short H2/H1 |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|\n');
for ic = 1:numel(results)
    for is = 1:numel(results(ic).fourier.station)
        s = results(ic).fourier.station(is);
        fprintf(fid, '| %s | %d | %.4g | %.4g | %.4g | %.4g |\n', ...
            results(ic).name, s.i, s.longAmpF0, s.shortAmpF0, s.longH2H1, s.shortH2H1);
    end
end

fprintf(fid, '\n## epsilon_fourier_summary.png\n\n');
fprintf(fid, '该图汇总平均真实主频幅值、归一化主频响应和二倍频比例。左图点旁标注的是未归一化的平均 `|P_f0|`，可直接看到 `E-2`、`E-3`、`E-4` 的真实压力扰动量级随 `epsilon` 逐级降低。中图除以 `epsilon` 后三组仍在相近量级，说明主频响应保留了线性缩放基础。右图的二倍频比例对非线性更敏感，`E-2` 明显高于小幅值算例，表明大扰动下波形畸变和高阶谐波生成增强。\n\n');

fprintf(fid, '## fourier_station_spectra_by_axis.png\n\n');
fprintf(fid, '该图为多站位 Fourier 频谱网格：三行分别对应 `i=%s`，左列为长轴，右列为短轴。每一行的长短轴子图使用相同纵轴范围，因此可直接比较同一站位下长轴和短轴谱峰高度。曲线按 `epsilon` 归一化，以便比较谱形；图例同时标注每条曲线的真实 `|P_f0|`，用于判断实际幅值大小。三组算例在 `f0` 附近均有主峰，说明后半段响应受外加周期扰动控制。沿流向看，后部站位的真实幅值通常减小；长短轴比较看，长轴主频幅值整体更大。`E-2` 在 `2f0` 附近的相对能量更强，说明大幅值算例包含更多非正弦成分。\n\n', join_ints(stationList));

fprintf(fid, '## fourier_station_i81_spectra.png\n\n');
fprintf(fid, '这是上一版保留下来的单站位频谱图，只比较 `i=81` 长轴方向的归一化频谱。本轮新版分析已经由 `fourier_station_spectra_by_axis.png` 取代它：新版同时包含 `i=%s` 三个站位，并分开展示长轴和短轴；因此论文或报告中建议优先使用新版多站位图。旧图仍可作为中部站位长轴频谱的快速核对图。\n\n', join_ints(stationList));

fprintf(fid, '## fourier_station_true_f0_amplitudes.png\n\n');
fprintf(fid, '该图用柱状图直接给出不同站位、长短轴上的真实主频幅值 `|P_f0|`，长轴和短轴两个面板使用相同纵轴范围，柱顶数字就是未归一化的压力扰动幅值。`E-2` 的真实主频幅值最大，`E-3` 次之，`E-4` 最小；同一算例中，长轴站位通常高于短轴站位，说明椭圆截面对快声波响应有稳定的周向放大差异。沿流向的站位差异反映了扰动传播、相位演化和衰减过程。\n\n');

fprintf(fid, '## fourier_station_second_harmonic_ratios.png\n\n');
fprintf(fid, '该图比较不同站位的 `|P_2f0|/|P_f0|`，长轴和短轴两个面板使用相同纵轴范围，柱顶数字为二倍频相对比例，图内文字给出真实 `|P_2f0|` 的范围。该图比主频幅值更适合判断非线性：若只是线性缩放，二倍频比例不会随 `epsilon` 显著升高。`E-2` 在多个站位上显著高于 `E-3` 和 `E-4`，说明大扰动幅值下二阶谐波生成明显增强。\n\n');

fprintf(fid, '## fourier_axis_f0_profiles.png\n\n');
fprintf(fid, '该图分为真实幅值和归一化幅值两类面板：左列给出长轴/短轴真实 `|P_f0|`，右列给出 `|P_f0|/epsilon`。上下对应的长轴和短轴面板保持相同纵轴范围，因此曲线高度可以直接比较。真实幅值面板显示 `E-2` 的实际压力扰动量级远高于小幅值算例；归一化面板显示三组曲线大体接近但不完全重合。长轴曲线整体高于短轴，说明椭圆锥体几何对入射扰动有稳定的周向调制作用。大幅值算例在局部位置的归一化偏离说明响应系数本身也受扰动幅值影响。\n\n');

fprintf(fid, '## fourier_axis_second_harmonic_ratio.png\n\n');
fprintf(fid, '该图同时给出真实二倍频幅值 `|P_2f0|` 和相对比例 `|P_2f0|/|P_f0|`。上下对应的长轴和短轴面板保持相同纵轴范围，便于直接判断哪一侧二倍频更强。真实二倍频幅值随 `epsilon` 增大而增大，而相对比例揭示波形畸变程度。`E-2` 的二倍频比例在长轴和短轴上都显著更高，是本轮分析中最清楚的非线性证据；`E-3` 和 `E-4` 的比例较低，整体更接近线性周期响应。\n\n');

fprintf(fid, '## pod_cumulative_energy.png\n\n');
fprintf(fid, 'POD 能量谱反映后半段壁面母线压力扰动能否由少数相干结构表示。图例标注了每个算例的 `epsilon` 和平均真实长轴 `|P_f0|`。前三组算例前两阶 POD 模态已达到 95%% 能量，说明主要动态可近似看作低维周期振荡。若大幅值算例的高阶能量占比增加，则通常对应 Fourier 中观察到的谐波增强和波形畸变。\n\n');

fprintf(fid, '## pod_long_axis_mode_shapes.png\n\n');
fprintf(fid, '该图比较前三个 POD 模态在长轴正向上的空间形状，图例中保留真实主频幅值量级。第一模态通常对应主导周期响应的主要空间包络，后续模态多与相位偏移、二倍频成分或局部畸变有关。若 `E-2` 的二、三阶模态在局部位置更突出，说明大幅值扰动改变的不只是整体强度，还改变了沿流向的空间结构。\n\n');

fprintf(fid, '## pod_time_coefficients.png\n\n');
fprintf(fid, 'POD 时间系数展示主要空间模态在后半段内如何随时间振荡。曲线本身已归一化以突出波形差异，图例给出真实幅值量级。若系数接近单频周期变化，说明该模态受基频控制；若大幅值算例的系数出现尖峰、非对称或波形畸变，则与 Fourier 中较高二倍频比例相互印证。\n\n');

fprintf(fid, '## dmd_eigenvalues.png\n\n');
fprintf(fid, 'DMD 特征值靠近单位圆表示对应模态在后半段时间窗中近似中性稳定、周期保持较好。图例中标出每个算例的真实主频幅值量级，便于把模态稳定性与扰动大小联系起来。三组算例的主要周期模态应靠近单位圆；若大幅值算例出现更多偏离单位圆的模态，说明它包含更复杂的瞬态或非线性调制成分。\n\n');

fprintf(fid, '## dmd_frequency_growth.png\n\n');
fprintf(fid, '该图把 DMD 模态映射到频率和增长率平面，点大小表示模态初始幅值，图例标注真实 `|P_f0|` 量级。靠近 `f0` 的大点代表主导周期响应，靠近 `2f0` 或更高频率的点代表高阶谐波。`E-2` 若在高频区域出现更大的点，说明二倍频及高阶成分不仅在 Fourier 中存在，也构成了可识别的动态模态。\n\n');

fprintf(fid, '## dmd_reconstruction_error.png\n\n');
fprintf(fid, '该图给出有限秩 DMD 对后半段快照矩阵的相对重构误差。横坐标标签同时含有 `epsilon` 和真实平均长轴 `|P_f0|`，因此可直接比较不同实际扰动量级下 DMD 的表示能力。误差越小，说明当前 DMD 秩越能表示准周期演化；误差较大通常意味着时间窗较短、采样点较少、包含非平稳成分，或需要更高秩模态描述波形畸变。\n\n');

fprintf(fid, '## dmd_f0_mode_shapes_long_axis.png\n\n');
fprintf(fid, '该图展示每个算例中最接近 `f0` 的 DMD 模态在长轴方向上的空间形状，标题中标注了真实幅值量级、对应频率和增长率。它可看作 DMD 意义下的主频传播结构。不同扰动幅值之间的模态形状若保持接近，说明主频响应的空间结构有线性相似性；局部偏离则说明幅值增大后，主频结构也受到非线性修正。\n\n');

fprintf(fid, '## 综合判断\n\n');
[~, idxMax] = max([results.epsilon]);
rmax = results(idxMax);
fprintf(fid, '综合 Fourier、POD 和 DMD，主频响应随扰动幅值增大而增强，且长轴方向始终强于短轴方向。归一化主频幅值说明响应仍保留线性缩放基础；二倍频比例、POD 高阶模态贡献和 DMD 高频模态则共同表明，大幅值 `%s` 已经出现更明显的非线性波形畸变。小幅值 `E-4` 更适合作为线性参考，中等幅值 `E-3` 介于两者之间。\n', rmax.name);
end

function sync_compare_outputs(rootDir, outDir, results, f0, stationList)
compareDir = fullfile(rootDir, 'Compare');
if ~exist(compareDir, 'dir')
    mkdir(compareDir);
end
figureNames = { ...
    'fourier_station_spectra_by_axis.png', ...
    'fourier_station_true_f0_amplitudes.png', ...
    'fourier_station_second_harmonic_ratios.png', ...
    'pod_long_axis_mode_shapes.png', ...
    'dmd_f0_mode_shapes_long_axis.png'};
for i = 1:numel(figureNames)
    copyfile(fullfile(outDir, figureNames{i}), fullfile(compareDir, figureNames{i}), 'f');
end
write_compare_report(compareDir, results, f0, stationList);
end

function write_compare_report(compareDir, results, f0, stationList)
fn = fullfile(compareDir, 'compare_figure_analysis.md');
fid = fopen_utf8(fn);
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Compare 图像润色与分析说明\n\n');
fprintf(fid, '本目录中的图已重新按论文展示风格绘制：画布增大，坐标轴、标题、图例字体放大，线宽加粗，图例尽量移到图外或面板上方，避免遮挡主要曲线。所有长轴/短轴并列对比图均保持成对一致的纵轴范围，便于直接比较幅值差异。\n\n');
fprintf(fid, '分析对象为壁面 `j=1` 的压力扰动 `P_pert`。扰动幅值分别为 `E-2: epsilon=5e-2`、`E-3: epsilon=5e-3`、`E-4: epsilon=5e-4`。Fourier 目标基频为 `f0=%.6g Hz`，站位为 `i=%s`。\n\n', f0, join_ints(stationList));

fprintf(fid, '## 关键数值表\n\n');
fprintf(fid, '| case | station | long P_f0 | short P_f0 | long H2/H1 | short H2/H1 |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|\n');
for ic = 1:numel(results)
    for is = 1:numel(results(ic).fourier.station)
        s = results(ic).fourier.station(is);
        fprintf(fid, '| %s | %d | %.4g | %.4g | %.4g | %.4g |\n', ...
            results(ic).name, s.i, s.longAmpF0, s.shortAmpF0, s.longH2H1, s.shortH2H1);
    end
end

fprintf(fid, '\n## fourier_station_spectra_by_axis.png\n\n');
fprintf(fid, '该图展示三个流向站位 `i=%s` 上长轴与短轴的后半段 Fourier 频谱。每一行的长轴和短轴面板使用相同纵轴范围，因此谱峰高度可以直接比较。曲线按 `epsilon` 归一化，图例中同时标出真实 `|P_f0|`，这样既能看谱形是否满足线性缩放，也能看真实压力扰动量级。\n\n', join_ints(stationList));
fprintf(fid, '三组算例在 `f0` 附近均存在清晰主峰，说明外加周期快声波是后半段响应的主导频率。前部站位 `i=40` 的主峰幅值最大，向 `i=81` 和 `i=121` 后移后明显减小，反映扰动沿流向传播时的衰减与相位演化。长轴面板的主频峰整体高于短轴，说明椭圆截面对入射扰动存在稳定的周向放大差异。`E-2` 在 `2f0` 附近的相对峰值更突出，是大幅值扰动引起波形畸变的重要证据。\n\n');

fprintf(fid, '## fourier_station_true_f0_amplitudes.png\n\n');
fprintf(fid, '该图直接比较各站位的真实主频幅值 `|P_f0|`，柱顶数字为未归一化压力扰动幅值。长轴和短轴面板使用相同纵轴范围，所以可以直接看出长轴响应通常高于短轴响应。真实量级上，`E-2` 最大、`E-3` 次之、`E-4` 最小，基本随扰动幅值逐级降低。\n\n');
fprintf(fid, '站位差异也很明显：`i=40` 处长轴真实主频幅值约为 `E-2: %.4g`、`E-3: %.4g`、`E-4: %.4g`，远高于中后部站位。这说明当前壁面压力扰动的主频响应主要集中在前部区域。中后部站位仍保持长轴略强于短轴的趋势，但绝对幅值已明显下降。\n\n', ...
    results(1).fourier.station(1).longAmpF0, results(2).fourier.station(1).longAmpF0, results(3).fourier.station(1).longAmpF0);

fprintf(fid, '## fourier_station_second_harmonic_ratios.png\n\n');
fprintf(fid, '该图比较各站位的二倍频比例 `|P_2f0|/|P_f0|`，长轴和短轴面板同样保持相同纵轴范围。二倍频比例越高，说明压力时间历程越偏离纯正弦主频响应，因而更能反映非线性波形畸变。\n\n');
fprintf(fid, '`E-2` 在多个站位上的二倍频比例明显高于 `E-3` 和 `E-4`，尤其在后部站位短轴方向可达到约 `%.4g`。这说明大幅值扰动不只是把小幅值响应整体放大，而是增强了高阶谐波。`E-3` 的二倍频比例整体较低，仍接近线性周期响应；`E-4` 是小幅值参考，但由于时间窗较短，个别站位的比例会受频谱泄漏和采样长度影响。\n\n', results(1).fourier.station(end).shortH2H1);

fprintf(fid, '## pod_long_axis_mode_shapes.png\n\n');
fprintf(fid, '该图比较前三个 POD 模态在长轴正向上的空间形状，图例放到图外以避免遮挡曲线。POD 第一模态通常代表后半段周期响应的主要空间包络，第二、第三模态则多与相位偏移、二倍频成分或局部波形畸变相关。\n\n');
fprintf(fid, '三组算例的一阶模态沿流向具有相似的主要结构，说明主频响应的空间分布有一定线性相似性。高阶模态更能体现幅值效应：大幅值 `E-2` 若在二、三阶模态上出现更明显的局部峰谷，说明非线性不只改变整体幅值，也会改变沿流向的空间结构。结合 Fourier 二倍频比例，可将这些高阶 POD 成分理解为波形畸变和高阶谐波在空间上的投影。\n\n');

fprintf(fid, '## dmd_f0_mode_shapes_long_axis.png\n\n');
fprintf(fid, '该图展示每个算例中最接近 `f0` 的 DMD 模态在长轴方向上的空间形状。标题中标出了对应频率和增长率，曲线表示 DMD 意义下的主频传播结构。由于 DMD 将时间演化和空间结构联系起来，该图可用于判断不同扰动幅值下主频动态结构是否保持相似。\n\n');
fprintf(fid, '三组算例都能提取到接近基频的 DMD 模态，说明后半段响应具有明确的周期动力学。若大幅值 `E-2` 的模态形状与小幅值 `E-4` 在局部区域出现偏离，可解释为主频结构受到非线性修正；若整体形状仍相近，则说明尽管二倍频增强，主频传播骨架仍保留了线性响应的主要特征。\n\n');

fprintf(fid, '## 综合结论\n\n');
fprintf(fid, '润色后的 Compare 图组共同表明：主频响应随 `epsilon` 增大而增大，且长轴方向整体强于短轴方向，前部站位响应最强。`E-2` 的二倍频比例显著升高，是当前不同扰动幅值对比中最清楚的非线性证据。POD 与 DMD 进一步说明，主频结构仍具有低维相干特征，但大幅值扰动会增强高阶成分并可能改变局部空间模态。\n');
end

function write_report(outDir, results, f0, stationI)
fn = fullfile(outDir, 'figure_analysis_report.md');
fid = fopen_utf8(fn);
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# 不同扰动幅值的 Fourier/POD/DMD 对比分析\n\n');
fprintf(fid, '本次后处理比较 `E-2`、`E-3`、`E-4` 三组快声波扰动幅值算例。分析对象为壁面 `j=1` 上四条代表母线：长轴正负向 `k=41,121` 与短轴正负向 `k=1,81` 的压力扰动 `P_pert`。为了避免初始瞬态主导模态，Fourier、POD 与 DMD 均主要使用每个算例后半段时间窗。\n\n');
fprintf(fid, '采样间隔为 `Delta t_save = %.6g`，目标基频取 `f0 = %.6g Hz`，频谱代表站位为长轴 `i = %d`。\n\n', results(1).dtSave, f0, stationI);

fprintf(fid, '## 数据与主要指标\n\n');
fprintf(fid, '| case | epsilon | frames | late frames | long P_f0/epsilon | short P_f0/epsilon | long H2/H1 | short H2/H1 | POD E1 | POD modes 95%% | DMD error |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ic = 1:numel(results)
    r = results(ic);
    fprintf(fid, '| %s | %.3g | %d | %d | %.4g | %.4g | %.4g | %.4g | %.4g | %d | %.4g |\n', ...
        r.name, r.epsilon, r.nFrames, numel(r.lateIdx), ...
        r.fourier.meanLongAmp / r.epsilon, r.fourier.meanShortAmp / r.epsilon, ...
        r.fourier.longH2H1, r.fourier.shortH2H1, r.pod.energy(1), r.pod.n95, r.dmd.relError);
end

fprintf(fid, '\n## epsilon_fourier_summary.png\n\n');
fprintf(fid, '该图把主频幅值、主频归一化响应和二倍频比例放在同一张图中。主频幅值随 `epsilon` 增大而增大，说明三组算例首先都响应在入射快声波的基频附近。归一化后，长轴方向大约保持在同一量级，短轴方向略低，说明壁面压力主频响应仍有明显线性缩放基础。二倍频比例则对非线性更敏感：`E-2` 明显高于小幅值算例，表明大扰动下波形畸变和高阶谐波生成增强。\n\n');

fprintf(fid, '## fourier_station_i81_spectra.png\n\n');
fprintf(fid, '该图比较长轴 `i=%d` 站位的后半段归一化频谱。三组算例在 `f0` 附近都有主峰，说明后半段已形成由外加周期扰动控制的稳定频率响应。大幅值 `E-2` 在二倍频及更高频率附近的能量相对更强，说明它不是简单地把小幅值结果乘以常数，而是产生了更多非正弦成分。\n\n', stationI);

fprintf(fid, '## fourier_axis_f0_profiles.png\n\n');
fprintf(fid, '该图给出长轴和短轴方向的 `|P_f0|/epsilon` 沿流向分布。长轴曲线整体高于短轴，说明椭圆锥体几何对入射扰动有稳定的周向调制作用。三组算例归一化后并不完全重合，尤其大幅值算例在局部位置会偏离小幅值参考，说明响应系数本身也受扰动幅值影响。\n\n');

fprintf(fid, '## fourier_axis_second_harmonic_ratio.png\n\n');
fprintf(fid, '该图沿流向展示 `|P_2f0|/|P_f0|`。二倍频比例越高，代表周期波形越偏离单纯正弦响应。`E-2` 的二倍频比例在长轴和短轴上都显著更高，是本轮分析中最清楚的非线性证据；`E-3` 和 `E-4` 的比例较低，整体更接近线性周期响应。\n\n');

fprintf(fid, '## pod_cumulative_energy.png\n\n');
fprintf(fid, 'POD 能量谱反映后半段壁面母线压力扰动能否由少数相干结构表示。若前一两个模态能量占比高，说明主要动态近似为低维周期振荡；如果需要更多模态才能达到 90%% 或 95%% 能量，则说明波形、相位或空间结构更复杂。大幅值算例通常会因谐波增强和波形畸变而提高高阶 POD 模态贡献。\n\n');

fprintf(fid, '## pod_long_axis_mode_shapes.png\n\n');
fprintf(fid, '该图比较前三个 POD 模态在长轴正向上的空间形状。第一模态通常对应主导周期响应的主要空间包络，后续模态多与相位偏移、二倍频成分或局部畸变有关。若 `E-2` 的二、三阶模态在局部位置更突出，说明大幅值扰动改变的不只是整体强度，还改变了沿流向的空间结构。\n\n');

fprintf(fid, '## pod_time_coefficients.png\n\n');
fprintf(fid, 'POD 时间系数展示主要空间模态在后半段内如何随时间振荡。若系数接近单频周期变化，说明该模态受基频控制；若大幅值算例的系数出现尖峰、非对称或波形畸变，则与 Fourier 中较高二倍频比例相互印证。\n\n');

fprintf(fid, '## dmd_eigenvalues.png\n\n');
fprintf(fid, 'DMD 特征值靠近单位圆表示对应模态在后半段时间窗中近似中性稳定、周期保持较好；偏离单位圆则对应增长或衰减。三组算例的主要周期模态应靠近单位圆，这与后半段准周期响应一致。若大幅值算例出现更多偏离单位圆的模态，说明它包含更复杂的瞬态或非线性调制成分。\n\n');

fprintf(fid, '## dmd_frequency_growth.png\n\n');
fprintf(fid, '该图把 DMD 模态映射到频率和增长率平面，点大小表示模态初始幅值。靠近 `f0` 的大幅值模态代表主导周期响应，靠近 `2f0` 或更高频率的模态代表高阶谐波。`E-2` 若在高频区域出现更大的点，说明二倍频及高阶成分不仅在 Fourier 中存在，也构成了可识别的动态模态。\n\n');

fprintf(fid, '## dmd_reconstruction_error.png\n\n');
fprintf(fid, '该图给出有限秩 DMD 对后半段快照矩阵的相对重构误差。误差越小，说明当前选取的 DMD 秩已经能较好表示准周期演化。误差较大通常意味着时间窗较短、采样点较少、包含非平稳成分，或者需要更高秩模态才能描述波形畸变。\n\n');

fprintf(fid, '## dmd_f0_mode_shapes_long_axis.png\n\n');
fprintf(fid, '该图展示每个算例中最接近 `f0` 的 DMD 模态在长轴方向上的空间形状。它可以看作 DMD 意义下的主频传播结构。不同扰动幅值之间的模态形状若保持接近，说明主频响应的空间结构有线性相似性；局部偏离则说明幅值增大后，主频结构也受到非线性修正。\n\n');

fprintf(fid, '## 综合判断\n\n');
[~, idxMax] = max([results.epsilon]);
rmax = results(idxMax);
fprintf(fid, '综合 Fourier、POD 和 DMD，主频响应随扰动幅值增大而增强，且长轴方向始终强于短轴方向。归一化主频幅值说明响应仍保留线性缩放基础；二倍频比例、POD 高阶模态贡献和 DMD 高频模态则共同表明，大幅值 `%s` 已经出现更明显的非线性波形畸变。小幅值 `E-4` 更适合作为线性参考，中等幅值 `E-3` 介于两者之间。\n', rmax.name);
end

function axesFrame = extract_wall_axes(masterPath, targetKs, needPoints)
pertDir = fileparts(masterPath);
info = parse_pvts(masterPath);
nAxes = numel(targetKs);
rows = cell(nAxes, 1);

for ip = 1:numel(info.pieces)
    ext = sscanf(info.pieces(ip).extent, '%f');
    if ext(3) > 1 || ext(4) < 1
        continue;
    end
    hitAxes = find(targetKs >= ext(5) & targetKs <= ext(6));
    if isempty(hitAxes)
        continue;
    end
    piecePath = fullfile(pertDir, info.pieces(ip).source);
    piece = read_vts_piece(piecePath, 'P_pert', needPoints);
    for ih = 1:numel(hitAxes)
        ia = hitAxes(ih);
        [idx, globI] = line_indices(piece.extent, 1, targetKs(ia));
        if needPoints
            x = piece.points(idx, 1);
        else
            x = nan(numel(idx), 1);
        end
        rows{ia} = [rows{ia}; globI(:), x(:), piece.values(idx)]; %#ok<AGROW>
    end
end

axesFrame = struct('i', {}, 'x', {}, 'p', {});
for ia = 1:nAxes
    if isempty(rows{ia})
        error('Cannot extract wall axis k=%d from %s.', targetKs(ia), masterPath);
    end
    [uI, ~, ic] = unique(rows{ia}(:, 1));
    pMean = accumarray(ic, rows{ia}(:, 3), [], @mean);
    if needPoints
        xMean = accumarray(ic, rows{ia}(:, 2), [], @mean);
    else
        xMean = nan(size(uI));
    end
    axesFrame(ia).i = uI; %#ok<AGROW>
    axesFrame(ia).x = xMean; %#ok<AGROW>
    axesFrame(ia).p = pMean; %#ok<AGROW>
end
end

function info = parse_pvts(pvtsFile)
txt = fileread(pvtsFile);
tok = regexp(txt, 'WholeExtent="\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)', 'tokens', 'once');
info.whole = str2double(tok);
ptoks = regexp(txt, '<Piece\s+Extent=\s*"\s*([^"]+)"\s+Source="([^"]+)"', 'tokens');
pieces = repmat(struct('extent', '', 'source', ''), 1, numel(ptoks));
for n = 1:numel(ptoks)
    pieces(n).extent = strtrim(ptoks{n}{1});
    pieces(n).source = ptoks{n}{2};
end
info.pieces = pieces;
end

function piece = read_vts_piece(pathName, varName, needPoints)
fid = fopen(pathName, 'r');
if fid < 0
    error('Cannot open %s.', pathName);
end
raw = fread(fid, '*uint8').';
fclose(fid);

txtAll = char(raw);
marker = '<AppendedData encoding="raw">_';
markerPos = strfind(txtAll, marker);
if isempty(markerPos)
    txt = txtAll;
else
    txt = txtAll(1:markerPos(1)-1);
end

extentTok = regexp(txt, '<Piece\s+Extent\s*=\s*"([^"]+)"', 'tokens', 'once');
if isempty(extentTok)
    error('Cannot find piece extent in %s.', pathName);
end
pieceExtent = extentTok{1};
nPoints = count_points_from_extent(pieceExtent);

points = [];
if needPoints
    pointsBlock = regexp(txt, '<Points>\s*(.*?)\s*</Points>', 'tokens', 'once');
    if isempty(pointsBlock)
        error('Cannot find Points block in %s.', pathName);
    end
    pointsArray = regexp(pointsBlock{1}, '<DataArray([^>]*)/?>', 'tokens', 'once');
    points = read_data_array(raw, txt, markerPos, marker, pointsArray{1}, pathName, 3*nPoints);
    points = reshape(points, 3, []).';
end

valueExpr = ['<DataArray([^>]*)Name="' regexptranslate('escape', varName) '"([^>]*)/?>'];
valueParts = regexp(txt, valueExpr, 'tokens', 'once');
if isempty(valueParts)
    error('Cannot find variable %s in %s.', varName, pathName);
end
values = read_data_array(raw, txt, markerPos, marker, [valueParts{1} ' ' valueParts{2}], pathName, nPoints);

piece.extent = pieceExtent;
piece.points = points;
piece.values = values;
end

function values = read_data_array(raw, headerText, markerPos, marker, attrs, pathName, expectedValues)
formatTok = regexp(attrs, 'format\s*=\s*"([^"]+)"', 'tokens', 'once');
format = 'ascii';
if ~isempty(formatTok)
    format = lower(strtrim(formatTok{1}));
end

if strcmp(format, 'appended')
    if isempty(markerPos)
        error('Appended data marker not found in %s.', pathName);
    end
    offsetTok = regexp(attrs, 'offset\s*=\s*"([^"]+)"', 'tokens', 'once');
    if isempty(offsetTok)
        error('Appended DataArray offset not found in %s.', pathName);
    end
    offset = str2double(strtrim(offsetTok{1}));
    dataStart = markerPos(1) + length(marker);
    blockStart = dataStart + offset;
    nBytes = expectedValues * 4;
    block = raw(blockStart+8:blockStart+8+nBytes-1);
    values = double(typecast(uint8(block), 'single')).';
else
    dataText = regexp(headerText, ['<DataArray' regexptranslate('escape', attrs) '[^>]*>(.*?)</DataArray>'], 'tokens', 'once');
    if isempty(dataText)
        error('ASCII DataArray contents not found in %s.', pathName);
    end
    values = sscanf(dataText{1}, '%f');
end
end

function [idx, globI] = line_indices(extentText, wallJ, targetK)
e = sscanf(extentText, '%f');
iVals = e(1):e(2);
jMin = e(3); jMax = e(4);
kMin = e(5); kMax = e(6);
if wallJ < jMin || wallJ > jMax || targetK < kMin || targetK > kMax
    idx = [];
    globI = [];
    return;
end
nI = e(2) - e(1) + 1;
nJ = e(4) - e(3) + 1;
jLocal = wallJ - jMin + 1;
kLocal = targetK - kMin + 1;
idx = (1:nI)' + (jLocal-1)*nI + (kLocal-1)*nI*nJ;
globI = iVals(:);
end

function nPoints = count_points_from_extent(extentText)
nums = sscanf(extentText, '%f');
nPoints = prod(nums(2:2:6) - nums(1:2:5) + 1);
end

function frameIds = extract_frame_ids(names)
frameIds = nan(numel(names), 1);
for i = 1:numel(names)
    tok = regexp(names{i}, '(\d+)\.pvts$', 'tokens', 'once');
    frameIds(i) = str2double(tok{1});
end
end

function cfg = read_config(fn)
txt = fileread(fn);
cfg.text = txt;
end

function val = get_cfg_number(cfg, key, fallback)
expr = ['(^|\n)\s*' key '\s*=\s*([+-]?\d*\.?\d+(?:[dDeE][+-]?\d+)?)'];
tok = regexp(cfg.text, expr, 'tokens', 'once');
if isempty(tok)
    val = fallback;
else
    val = str2double(regexprep(tok{2}, '[dD]', 'e'));
    if ~isfinite(val)
        val = fallback;
    end
end
end

function [X, meta] = stack_axes(pAxes, iAxes, xAxes, axisNames)
rows = [];
axisId = [];
iId = [];
xVal = [];
axisLabel = {};
for ia = 1:numel(pAxes)
    rows = [rows; pAxes{ia}]; %#ok<AGROW>
    n = size(pAxes{ia}, 1);
    axisId = [axisId; ia * ones(n, 1)]; %#ok<AGROW>
    iId = [iId; iAxes{ia}(:)]; %#ok<AGROW>
    xVal = [xVal; xAxes{ia}(:)]; %#ok<AGROW>
    axisLabel = [axisLabel; repmat(axisNames(ia), n, 1)]; %#ok<AGROW>
end
X = rows;
meta.axisId = axisId;
meta.i = iId;
meta.x = xVal;
meta.axisLabel = axisLabel;
end

function X = fill_nan_matrix(X)
for i = 1:size(X, 1)
    row = X(i, :);
    ok = isfinite(row);
    if all(ok)
        continue;
    elseif any(ok)
        row(~ok) = mean(row(ok));
    else
        row(:) = 0;
    end
    X(i, :) = row;
end
end

function a = mean_same_grid(cells)
a = zeros(size(cells{1}));
count = zeros(size(cells{1}));
for i = 1:numel(cells)
    v = cells{i};
    ok = isfinite(v);
    a(ok) = a(ok) + v(ok);
    count(ok) = count(ok) + 1;
end
a = a ./ max(count, 1);
end

function m = circular_mean(cells)
z = zeros(size(cells{1}));
count = zeros(size(cells{1}));
for i = 1:numel(cells)
    v = cells{i};
    ok = isfinite(v);
    z(ok) = z(ok) + exp(1i * v(ok));
    count(ok) = count(ok) + 1;
end
m = unwrap(angle(z ./ max(count, 1)));
end

function [freq, amp, raw] = station_spectrum(pAxes, iAxes, lateIdx, t, fs, targetI, axisRange)
sig = [];
for ia = axisRange
    [~, loc] = min(abs(iAxes{ia} - targetI));
    sig = [sig; pAxes{ia}(loc, :)]; %#ok<AGROW>
end
y = mean(sig, 1, 'omitnan');
y = y(lateIdx);
y = y - mean(y, 'omitnan');
y = y(:);
n = numel(y);
if n < 4
    freq = 0;
    amp = 0;
    raw = y;
    return;
end
w = local_hann(n);
nfft = 4096;
Y = fft(y .* w, nfft);
nHalf = floor(nfft/2) + 1;
freq = (0:nHalf-1)' * fs / nfft;
amp = 2 * abs(Y(1:nHalf)) / sum(w);
amp(1) = 0;
raw = y;
end

function station = station_fourier_metrics(pAxes, iAxes, lateIdx, t, fs, stationList, f0)
station = repmat(struct( ...
    'i', NaN, ...
    'longAmpF0', NaN, 'shortAmpF0', NaN, ...
    'longAmp2F0', NaN, 'shortAmp2F0', NaN, ...
    'longH2H1', NaN, 'shortH2H1', NaN, ...
    'freqLong', [], 'specLong', [], ...
    'freqShort', [], 'specShort', []), numel(stationList), 1);
for is = 1:numel(stationList)
    targetI = stationList(is);
    [freqLong, specLong, rawLong] = station_spectrum(pAxes, iAxes, lateIdx, t, fs, targetI, 1:2);
    [freqShort, specShort, rawShort] = station_spectrum(pAxes, iAxes, lateIdx, t, fs, targetI, 3:4);
    station(is).i = targetI;
    station(is).freqLong = freqLong;
    station(is).specLong = specLong;
    station(is).freqShort = freqShort;
    station(is).specShort = specShort;
    station(is).longAmpF0 = harmonic_amplitude(rawLong, t(lateIdx), f0);
    station(is).shortAmpF0 = harmonic_amplitude(rawShort, t(lateIdx), f0);
    station(is).longAmp2F0 = harmonic_amplitude(rawLong, t(lateIdx), 2*f0);
    station(is).shortAmp2F0 = harmonic_amplitude(rawShort, t(lateIdx), 2*f0);
    station(is).longH2H1 = station(is).longAmp2F0 / max(station(is).longAmpF0, eps);
    station(is).shortH2H1 = station(is).shortAmp2F0 / max(station(is).shortAmpF0, eps);
end
end

function amp = harmonic_amplitude(y, t, freq)
y = y(:);
t = t(:);
ok = isfinite(y) & isfinite(t);
y = y(ok);
t = t(ok);
if numel(y) < 2
    amp = NaN;
    return;
end
y = y - mean(y);
coeff = (2 / numel(y)) * (y.' * exp(-1i * 2*pi*freq*t));
amp = abs(coeff);
end

function w = local_hann(n)
if n <= 1
    w = ones(n, 1);
else
    i = (0:n-1)';
    w = 0.5 - 0.5*cos(2*pi*i/(n-1));
end
end

function [x, y] = axis_component(vec, meta, axisNumber)
idx = meta.axisId == axisNumber;
x = meta.x(idx);
y = vec(idx);
[x, order] = sort(x);
y = y(order);
end

function labels = case_labels(results)
labels = cell(1, numel(results));
for i = 1:numel(results)
    labels{i} = sprintf('%s eps=%.1e, <|P_f0|>L=%.3g', ...
        results(i).name, results(i).epsilon, results(i).fourier.meanLongAmp);
end
end

function labels = case_names(results)
labels = cell(1, numel(results));
for i = 1:numel(results)
    labels{i} = sprintf('%s, eps=%.1e', results(i).name, results(i).epsilon);
end
end

function yMax = station_spectrum_ymax(results, stationIndex)
yMax = 0;
for ic = 1:numel(results)
    s = results(ic).fourier.station(stationIndex);
    yMax = max(yMax, max(s.specLong ./ results(ic).epsilon));
    yMax = max(yMax, max(s.specShort ./ results(ic).epsilon));
end
yMax = padded_ymax(yMax);
end

function yMax = axis_f0_ymax(results, normalized)
yMax = 0;
for ic = 1:numel(results)
    if normalized
        longVals = results(ic).fourier.longAmp ./ results(ic).epsilon;
        shortVals = results(ic).fourier.shortAmp ./ results(ic).epsilon;
    else
        longVals = results(ic).fourier.longAmp;
        shortVals = results(ic).fourier.shortAmp;
    end
    yMax = max(yMax, max(longVals));
    yMax = max(yMax, max(shortVals));
end
yMax = padded_ymax(yMax);
end

function yMax = axis_h2_ymax(results, ratioMode)
yMax = 0;
for ic = 1:numel(results)
    if ratioMode
        longVals = results(ic).fourier.longAmp2 ./ max(results(ic).fourier.longAmp, eps);
        shortVals = results(ic).fourier.shortAmp2 ./ max(results(ic).fourier.shortAmp, eps);
    else
        longVals = results(ic).fourier.longAmp2;
        shortVals = results(ic).fourier.shortAmp2;
    end
    yMax = max(yMax, max(longVals));
    yMax = max(yMax, max(shortVals));
end
yMax = padded_ymax(yMax);
end

function yMax = padded_ymax(yMax)
if ~isfinite(yMax) || yMax <= 0
    yMax = 1;
else
    yMax = yMax * 1.12;
end
end

function txt = join_ints(vals)
parts = cell(1, numel(vals));
for i = 1:numel(vals)
    parts{i} = sprintf('%d', vals(i));
end
txt = strjoin(parts, ', ');
end

function annotate_points(x, y, values)
for i = 1:numel(x)
    text(x(i), y(i), sprintf(' %.3g', values(i)), ...
        'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end
end

function annotate_grouped_bars(vals, yMax)
[nGroups, nSeries] = size(vals);
groupWidth = min(0.8, nSeries/(nSeries + 1.5));
yMax = padded_ymax(yMax);
for i = 1:nGroups
    for j = 1:nSeries
        x = i - groupWidth/2 + (2*j-1) * groupWidth / (2*nSeries);
        text(x, vals(i, j) + 0.025*yMax, sprintf('%.3g', vals(i, j)), ...
            'Rotation', 90, 'FontSize', 9, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end
ylim([0, yMax]);
end

function v = get_or_nan(x, i)
if numel(x) >= i
    v = x(i);
else
    v = NaN;
end
end

function fig = make_figure(pos)
fig = figure('Visible', 'off', 'Position', [100, 100, pos], 'Color', 'w');
end

function style_pub(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 22, 'FontWeight', 'bold', 'LineWidth', 1.6);
set(get(ax, 'XLabel'), 'FontName', 'Times New Roman', 'FontSize', 25, 'FontWeight', 'bold');
set(get(ax, 'YLabel'), 'FontName', 'Times New Roman', 'FontSize', 25, 'FontWeight', 'bold');
set(get(ax, 'Title'), 'FontName', 'Times New Roman', 'FontSize', 27, 'FontWeight', 'bold');
box(ax, 'on');
end

function save_fig(fig, fn)
print(fig, fn, '-dpng', '-r300');
close(fig);
end

function fid = fopen_utf8(fn)
fid = fopen(fn, 'w', 'n', 'UTF-8');
if fid < 0
    fid = fopen(fn, 'w');
end
end
