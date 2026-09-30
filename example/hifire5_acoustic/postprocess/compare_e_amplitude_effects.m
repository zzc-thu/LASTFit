function compare_e_amplitude_effects()
% Compare E-2/E-3/E-4 fast-wave perturbation amplitudes using Pert fields.

scriptDir = fileparts(mfilename('fullpath'));
caseReleaseDir = fileparts(scriptDir);
rootDir = getenv('LASTFIT_HIFIRE_RAW_DIR');
if isempty(rootDir)
    error(['Set LASTFIT_HIFIRE_RAW_DIR to the directory containing ', ...
        'E-2, E-3, E-4, and RESU (for example, Un_Hifire_fine).']);
end
cases = struct( ...
    'name', {'E-2', 'E-3', 'E-4'}, ...
    'epsilon', {5e-2, 5e-3, 5e-4}, ...
    'epsilonSource', {'inferred from folder name', 'inferred from folder name', 'Config.cfg'});

outDir = fullfile(caseReleaseDir, 'output', 'manuscript_fine', ...
    'amplitude_comparison');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

dtSave = 2e-4;
fs = 1 / dtSave;
f0 = 666.6666666667;
targetKs = [41, 121, 1, 81];
axisNames = {'Long axis +', 'Long axis -', 'Short axis +', 'Short axis -'};
stationI = [40, 81, 121];

results = struct([]);
for ic = 1:numel(cases)
    fprintf('\n===== Processing %s =====\n', cases(ic).name);
    caseDir = fullfile(rootDir, cases(ic).name);
    cfgPath = fullfile(caseDir, 'Config.cfg');
    if exist(cfgPath, 'file')
        cfg = read_config(cfgPath);
        cases(ic).epsilon = get_cfg_number(cfg, 'epsilon', cases(ic).epsilon);
        cases(ic).epsilonSource = 'Config.cfg';
        dt = get_cfg_number(cfg, 'FixedDeltaTime', 1e-7);
        steps = get_cfg_number(cfg, 'StepsWriteData', 2000);
        dtSave = dt * steps;
        fs = 1 / dtSave;
    end

    result = analyze_case(caseDir, cases(ic), dtSave, fs, f0, targetKs, axisNames, stationI);
    results = [results; result]; %#ok<AGROW>
end

write_summary(outDir, results, f0, dtSave, fs);
write_metrics_csv(outDir, results);
plot_scaling(outDir, results);
plot_normalized_response(outDir, results);
plot_second_harmonic(outDir, results);
plot_axis_profiles(outDir, results);
plot_phase_profiles(outDir, results);
plot_time_histories(outDir, results, stationI);
plot_spectra(outDir, results, f0);
plot_spectra_raw(outDir, results, f0);
plot_nonlinearity_index(outDir, results);

fprintf('\nDone. Output directory:\n%s\n', outDir);
end

function result = analyze_case(caseDir, caseInfo, dtSave, fs, f0, targetKs, axisNames, stationI)
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
t = (0:nFrames-1) * dtSave;
lateIdx = ceil(nFrames/2):nFrames;
tLate = t(lateIdx);

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

f2 = 2 * f0;
nLate = numel(lateIdx);
phase1 = exp(-1i * 2*pi*f0*tLate(:));
phase2 = exp(-1i * 2*pi*f2*tLate(:));

amp1 = cell(nAxes, 1);
amp2 = cell(nAxes, 1);
phaseRad = cell(nAxes, 1);
rmsLate = cell(nAxes, 1);
for ia = 1:nAxes
    y = pAxes{ia}(:, lateIdx);
    y = y - mean(y, 2, 'omitnan');
    coeff1 = (2 / nLate) * (y * phase1);
    coeff2 = (2 / nLate) * (y * phase2);
    amp1{ia} = abs(coeff1);
    amp2{ia} = abs(coeff2);
    phaseRad{ia} = unwrap(angle(coeff1));
    rmsLate{ia} = sqrt(mean(y.^2, 2, 'omitnan'));
end

longIdx = 1:2;
shortIdx = 3:4;
xCommon = xAxes{1};
iCommon = iAxes{1};
validMask = iCommon >= 10;

longAmp = mean_same_grid(amp1(longIdx));
shortAmp = mean_same_grid(amp1(shortIdx));
longAmp2 = mean_same_grid(amp2(longIdx));
shortAmp2 = mean_same_grid(amp2(shortIdx));
longRms = mean_same_grid(rmsLate(longIdx));
shortRms = mean_same_grid(rmsLate(shortIdx));
longPhase = circular_mean(phaseRad(longIdx));
shortPhase = circular_mean(phaseRad(shortIdx));

meanLongAmp = mean(longAmp(validMask), 'omitnan');
meanShortAmp = mean(shortAmp(validMask), 'omitnan');
meanLongAmp2 = mean(longAmp2(validMask), 'omitnan');
meanShortAmp2 = mean(shortAmp2(validMask), 'omitnan');
meanLongRms = mean(longRms(validMask), 'omitnan');
meanShortRms = mean(shortRms(validMask), 'omitnan');

[dominantFrequency, spectrum] = station_spectrum(pAxes, iAxes, t, fs, stationI(2));
station = station_metrics(pAxes, iAxes, lateIdx, tLate, f0, stationI);

result.name = caseInfo.name;
result.epsilon = caseInfo.epsilon;
result.epsilonSource = caseInfo.epsilonSource;
result.caseDir = caseDir;
result.frameIds = frameIds;
result.nFrames = nFrames;
result.lateIdx = lateIdx;
result.t = t;
result.tLate = tLate;
result.dtSave = dtSave;
result.fs = fs;
result.targetKs = targetKs;
result.axisNames = axisNames;
result.i = iCommon;
result.x = xCommon;
result.pAxes = pAxes;
result.iAxes = iAxes;
result.xAxes = xAxes;
result.amp1 = amp1;
result.amp2 = amp2;
result.phaseRad = phaseRad;
result.rmsLate = rmsLate;
result.longAmp = longAmp;
result.shortAmp = shortAmp;
result.longAmp2 = longAmp2;
result.shortAmp2 = shortAmp2;
result.longRms = longRms;
result.shortRms = shortRms;
result.longPhase = longPhase;
result.shortPhase = shortPhase;
result.meanLongAmp = meanLongAmp;
result.meanShortAmp = meanShortAmp;
result.meanLongAmp2 = meanLongAmp2;
result.meanShortAmp2 = meanShortAmp2;
result.meanLongRms = meanLongRms;
result.meanShortRms = meanShortRms;
result.meanLongH2H1 = meanLongAmp2 / meanLongAmp;
result.meanShortH2H1 = meanShortAmp2 / meanShortAmp;
result.longShortAmpRatio = meanLongAmp / meanShortAmp;
result.meanLongRms = meanLongRms;
result.meanShortRms = meanShortRms;
result.dominantFrequency = dominantFrequency;
result.spectrum = spectrum;
result.station = station;
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

function nPoints = count_points_from_extent(extentText)
nums = sscanf(extentText, '%f');
nPoints = prod(nums(2:2:6) - nums(1:2:5) + 1);
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

function info = parse_pvts(pathName)
txt = fileread(pathName);
whole = regexp(txt, 'WholeExtent="([^"]+)"', 'tokens', 'once');
if isempty(whole)
    error('Cannot find WholeExtent in %s.', pathName);
end
pieceTokens = regexp(txt, '<Piece\s+Extent\s*=\s*"([^"]+)"\s+Source\s*=\s*"([^"]+)"\s*/?>', 'tokens');
if isempty(pieceTokens)
    error('Cannot find Piece entries in %s.', pathName);
end
info.wholeExtent = whole{1};
info.pieces = struct('extent', {}, 'source', {});
for i = 1:numel(pieceTokens)
    info.pieces(i).extent = pieceTokens{i}{1}; %#ok<AGROW>
    info.pieces(i).source = pieceTokens{i}{2}; %#ok<AGROW>
end
end

function ids = extract_frame_ids(names)
ids = zeros(numel(names), 1);
for i = 1:numel(names)
    tok = regexp(names{i}, '(?:Pert_)?Result_(\d+)\.pvts', 'tokens', 'once');
    if isempty(tok)
        ids(i) = i;
    else
        ids(i) = str2double(tok{1});
    end
end
end

function cfg = read_config(pathName)
cfg = containers.Map();
txt = fileread(pathName);
lines = regexp(txt, '\r?\n', 'split');
for i = 1:numel(lines)
    line = strtrim(regexprep(lines{i}, '[#!%].*$', ''));
    tok = regexp(line, '^\s*([^=\s]+)\s*=\s*([^,\s]+)', 'tokens', 'once');
    if ~isempty(tok)
        cfg(tok{1}) = tok{2};
    end
end
end

function value = get_cfg_number(cfg, key, defaultValue)
value = defaultValue;
if isKey(cfg, key)
    tmp = str2double(cfg(key));
    if ~isnan(tmp)
        value = tmp;
    end
end
end

function out = mean_same_grid(cellsIn)
mat = nan(numel(cellsIn{1}), numel(cellsIn));
for i = 1:numel(cellsIn)
    mat(:, i) = cellsIn{i};
end
out = mean(mat, 2, 'omitnan');
end

function out = circular_mean(phaseCells)
mat = nan(numel(phaseCells{1}), numel(phaseCells));
for i = 1:numel(phaseCells)
    mat(:, i) = phaseCells{i};
end
out = angle(mean(exp(1i*mat), 2, 'omitnan'));
out = unwrap(out);
end

function [dominantFrequency, spectrum] = station_spectrum(pAxes, iAxes, t, fs, targetI)
sig = station_signal(pAxes, iAxes, targetI, 1:2);
sig = mean(sig, 2, 'omitnan');
sig = sig - mean(sig, 'omitnan');
n = numel(sig);
y = fft(sig);
nHalf = floor(n/2) + 1;
freq = (0:nHalf-1)' * fs / n;
amp = 2 * abs(y(1:nHalf)) / n;
amp(1) = 0;
[~, imax] = max(amp);
dominantFrequency = freq(imax);
spectrum.freq = freq;
spectrum.amp = amp;
spectrum.signal = sig;
spectrum.time = t(:);
end

function station = station_metrics(pAxes, iAxes, lateIdx, tLate, f0, stationI)
station = struct([]);
for is = 1:numel(stationI)
    longSig = mean(station_signal(pAxes, iAxes, stationI(is), 1:2), 2, 'omitnan');
    shortSig = mean(station_signal(pAxes, iAxes, stationI(is), 3:4), 2, 'omitnan');
    longLate = longSig(lateIdx);
    shortLate = shortSig(lateIdx);
    longLate = longLate - mean(longLate, 'omitnan');
    shortLate = shortLate - mean(shortLate, 'omitnan');
    cLong = 2/numel(tLate) * (longLate(:).' * exp(-1i*2*pi*f0*tLate(:)));
    cShort = 2/numel(tLate) * (shortLate(:).' * exp(-1i*2*pi*f0*tLate(:)));
    station(is).i = stationI(is); %#ok<AGROW>
    station(is).longAmp = abs(cLong); %#ok<AGROW>
    station(is).shortAmp = abs(cShort); %#ok<AGROW>
    station(is).longPhase = angle(cLong); %#ok<AGROW>
    station(is).shortPhase = angle(cShort); %#ok<AGROW>
    station(is).longPeakToPeak = max(longLate) - min(longLate); %#ok<AGROW>
    station(is).shortPeakToPeak = max(shortLate) - min(shortLate); %#ok<AGROW>
end
end

function sig = station_signal(pAxes, iAxes, targetI, axisIdx)
sig = nan(size(pAxes{axisIdx(1)}, 2), numel(axisIdx));
for j = 1:numel(axisIdx)
    ia = axisIdx(j);
    [~, loc] = min(abs(iAxes{ia} - targetI));
    sig(:, j) = pAxes{ia}(loc, :).';
end
end

function plot_scaling(outDir, results)
fig = make_figure();
epsVals = [results.epsilon];
long = [results.meanLongAmp];
short = [results.meanShortAmp];
loglog(epsVals, long, '-o', 'LineWidth', 2.2, 'MarkerSize', 7); hold on;
loglog(epsVals, short, '-s', 'LineWidth', 2.2, 'MarkerSize', 7);
[epsSort, idx] = sort(epsVals);
ref = long(idx(1)) * epsSort / epsSort(1);
loglog(epsSort, ref, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('\epsilon');
ylabel('|P_{f0}|');
legend('Long axis', 'Short axis', 'Linear reference', 'Location', 'northwest');
title('Fundamental harmonic amplitude scaling');
save_fig(fig, fullfile(outDir, 'fundamental_amplitude_scaling.png'));
end

function plot_normalized_response(outDir, results)
fig = make_figure();
epsVals = [results.epsilon];
semilogx(epsVals, [results.meanLongAmp] ./ epsVals, '-o', 'LineWidth', 2.2, 'MarkerSize', 7); hold on;
semilogx(epsVals, [results.meanShortAmp] ./ epsVals, '-s', 'LineWidth', 2.2, 'MarkerSize', 7);
grid on;
xlabel('\epsilon');
ylabel('|P_{f0}| / \epsilon');
legend('Long axis', 'Short axis', 'Location', 'best');
title('Normalized fundamental response');
save_fig(fig, fullfile(outDir, 'normalized_fundamental_response.png'));
end

function plot_second_harmonic(outDir, results)
fig = make_figure();
epsVals = [results.epsilon];
semilogx(epsVals, [results.meanLongH2H1], '-o', 'LineWidth', 2.2, 'MarkerSize', 7); hold on;
semilogx(epsVals, [results.meanShortH2H1], '-s', 'LineWidth', 2.2, 'MarkerSize', 7);
grid on;
xlabel('\epsilon');
ylabel('|P_{2f0}| / |P_{f0}|');
legend('Long axis', 'Short axis', 'Location', 'best');
title('Second-harmonic content as a nonlinearity indicator');
save_fig(fig, fullfile(outDir, 'second_harmonic_ratio.png'));
end

function plot_axis_profiles(outDir, results)
fig = make_figure([1200 760]);
colors = lines(numel(results));
subplot(2,1,1);
for ic = 1:numel(results)
    plot(results(ic).x, results(ic).longAmp ./ results(ic).epsilon, '-', 'Color', colors(ic,:), 'LineWidth', 1.9); hold on;
    plot(results(ic).x, results(ic).shortAmp ./ results(ic).epsilon, '--', 'Color', colors(ic,:), 'LineWidth', 1.9);
end
grid on;
xlabel('X');
ylabel('|P_{f0}| / \epsilon');
title('Normalized f0 amplitude profile');
legend(case_legend(results), 'Location', 'best');
subplot(2,1,2);
for ic = 1:numel(results)
    plot(results(ic).x, results(ic).longAmp, '-', 'Color', colors(ic,:), 'LineWidth', 1.9); hold on;
    plot(results(ic).x, results(ic).shortAmp, '--', 'Color', colors(ic,:), 'LineWidth', 1.9);
end
grid on;
xlabel('X');
ylabel('|P_{f0}|');
title('Raw f0 amplitude profile');
save_fig(fig, fullfile(outDir, 'axis_f0_amplitude_profiles.png'));
end

function plot_phase_profiles(outDir, results)
fig = make_figure([1200 680]);
colors = lines(numel(results));
for ic = 1:numel(results)
    plot(results(ic).x, results(ic).longPhase, '-', 'Color', colors(ic,:), 'LineWidth', 1.9); hold on;
    plot(results(ic).x, results(ic).shortPhase, '--', 'Color', colors(ic,:), 'LineWidth', 1.9);
end
grid on;
xlabel('X');
ylabel('phase(P_{f0}) [rad]');
title('Fundamental phase profile');
legend(case_legend(results), 'Location', 'best');
save_fig(fig, fullfile(outDir, 'axis_f0_phase_profiles.png'));
end

function plot_time_histories(outDir, results, stationI)
fig = make_figure([1300 900]);
colors = lines(numel(results));
f0 = 666.6666666667;
displayEnd = 3 / f0;
for is = 1:numel(stationI)
    ax = subplot(numel(stationI), 1, is);
    for ic = 1:numel(results)
        sig = mean(station_signal(results(ic).pAxes, results(ic).iAxes, stationI(is), 1:2), 2, 'omitnan');
        late = results(ic).lateIdx;
        tt = results(ic).t(late);
        tt = tt - tt(1);
        yy = sig(late) / results(ic).epsilon;
        keep = tt <= displayEnd + results(ic).dtSave/2;
        ttPlot = tt(keep);
        yyPlot = yy(keep);
        [ts, ys] = periodic_fit_time_series(tt(:), yy(:), f0, displayEnd);
        plot(ts, ys, '-', 'Color', colors(ic,:), 'LineWidth', 2.2); hold on;
        plot(ttPlot, yyPlot, 'o', 'Color', colors(ic,:), 'MarkerSize', 4.2, ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    grid on;
    xlabel('t - t_{late,start}');
    ylabel('P'' / \epsilon');
    title(sprintf('Long-axis normalized late time history, i = %d', stationI(is)));
    xlim([0, displayEnd]);
    if is == 1
        lgd = legend({results.name}, 'Location', 'best');
        set(lgd, 'FontName','Times New Roman','FontSize',22,'FontWeight','bold');
    end
    set(ax, 'FontName','Times New Roman','FontSize',18,'FontWeight','bold','LineWidth',1.2);
end
save_fig(fig, fullfile(outDir, 'late_time_histories_normalized_long_axis.png'));
end

function [ts, ys] = periodic_fit_time_series(t, y, f0, displayEnd)
t = t(:);
y = y(:);
ok = isfinite(t) & isfinite(y);
t = t(ok);
y = y(ok);
ts = linspace(0, displayEnd, 900).';
if numel(t) < 6
    ys = interp1(t, y, ts, 'pchip', 'extrap');
    return;
end
omega = 2*pi*f0;
H = [ones(numel(t),1), cos(omega*t), sin(omega*t), ...
     cos(2*omega*t), sin(2*omega*t)];
coef = H \ y;
ys = [ones(numel(ts),1), cos(omega*ts), sin(omega*ts), ...
      cos(2*omega*ts), sin(2*omega*ts)] * coef;
ys = min(max(ys, min(y)), max(y));
end

function plot_spectra(outDir, results, f0)
fig = make_figure([1200 720]);
colors = lines(numel(results));
commonEnd = common_late_duration(results);
for ic = 1:numel(results)
    [freq, amp] = common_late_station_spectrum(results(ic), 81, commonEnd);
    plot(freq, amp ./ results(ic).epsilon, '-', ...
        'Color', colors(ic,:), 'LineWidth', 2.0); hold on;
end
xline(f0, 'k--', 'LineWidth', 1.6);
xline(2*f0, 'k:', 'LineWidth', 1.6);
grid on;
xlabel('frequency');
ylabel('spectrum amplitude / \epsilon');
title('Late-window long-axis pressure spectrum at i = 81');
lgd = legend({results.name}, 'Location', 'best');
set(lgd, 'FontName','Times New Roman','FontSize',22,'FontWeight','bold');
xlim([0, 2200]);
set(gca, 'FontName','Times New Roman','FontSize',18,'FontWeight','bold','LineWidth',1.2);
save_fig(fig, fullfile(outDir, 'station_i81_spectra_normalized.png'));
end

function plot_spectra_raw(outDir, results, f0)
fig = make_figure([1200 720]);
colors = lines(numel(results));
commonEnd = common_late_duration(results);
for ic = 1:numel(results)
    [freq, amp] = common_late_station_spectrum_raw(results(ic), 81, commonEnd);
    plot(freq, amp ./ results(ic).epsilon, '-o', ...
        'Color', colors(ic,:), 'LineWidth', 1.6, 'MarkerSize', 5.0); hold on;
end
xline(f0, 'k--', 'LineWidth', 1.6);
xline(2*f0, 'k:', 'LineWidth', 1.6);
grid on;
xlabel('frequency');
ylabel('spectrum amplitude / \epsilon');
title('Raw late-window long-axis pressure spectrum at i = 81');
lgd = legend({results.name}, 'Location', 'best');
set(lgd, 'FontName','Times New Roman','FontSize',22,'FontWeight','bold');
xlim([0, 2200]);
set(gca, 'FontName','Times New Roman','FontSize',18,'FontWeight','bold','LineWidth',1.2);
save_fig(fig, fullfile(outDir, 'station_i81_spectra_normalized_raw.png'));
end

function duration = common_late_duration(results)
duration = inf;
for ic = 1:numel(results)
    late = results(ic).lateIdx;
    tt = results(ic).t(late);
    tt = tt - tt(1);
    duration = min(duration, max(tt));
end
end

function [freq, amp] = common_late_station_spectrum(result, targetI, commonEnd)
sig = mean(station_signal(result.pAxes, result.iAxes, targetI, 1:2), 2, 'omitnan');
late = result.lateIdx;
tLate = result.t(late);
tRel = tLate - tLate(1);
keep = tRel <= commonEnd + 10*eps(commonEnd);
y = sig(late(keep));
y = y - mean(y, 'omitnan');
y = y(:);
n = numel(y);
if n < 4
    freq = 0;
    amp = 0;
    return;
end
w = local_hann(n);
nfft = 4096;
Y = fft(y .* w, nfft);
nHalf = floor(nfft/2) + 1;
freq = (0:nHalf-1)' * result.fs / nfft;
amp = 2 * abs(Y(1:nHalf)) / sum(w);
amp(1) = 0;
end

function [freq, amp] = common_late_station_spectrum_raw(result, targetI, commonEnd)
sig = mean(station_signal(result.pAxes, result.iAxes, targetI, 1:2), 2, 'omitnan');
late = result.lateIdx;
tLate = result.t(late);
tRel = tLate - tLate(1);
keep = tRel <= commonEnd + 10*eps(commonEnd);
y = sig(late(keep));
y = y - mean(y, 'omitnan');
y = y(:);
n = numel(y);
if n < 2
    freq = 0;
    amp = 0;
    return;
end
Y = fft(y);
nHalf = floor(n/2) + 1;
freq = (0:nHalf-1)' * result.fs / n;
amp = 2 * abs(Y(1:nHalf)) / n;
amp(1) = 0;
end

function w = local_hann(n)
if n <= 1
    w = ones(n,1);
else
    k = (0:n-1)';
    w = 0.5 - 0.5*cos(2*pi*k/(n-1));
end
end

function plot_nonlinearity_index(outDir, results)
epsVals = [results.epsilon];
[~, idxRef] = min(epsVals);
baseLong = results(idxRef).meanLongAmp / results(idxRef).epsilon;
baseShort = results(idxRef).meanShortAmp / results(idxRef).epsilon;
devLong = ([results.meanLongAmp] ./ epsVals) / baseLong - 1;
devShort = ([results.meanShortAmp] ./ epsVals) / baseShort - 1;

fig = make_figure([1100 620]);
bar(categorical({results.name}), [devLong(:), devShort(:)] * 100);
grid on;
ylabel('deviation from smallest-\epsilon linear reference [%]');
title('Fundamental-response nonlinearity index');
legend('Long axis', 'Short axis', 'Location', 'best');
save_fig(fig, fullfile(outDir, 'fundamental_nonlinearity_index.png'));
end

function labels = case_legend(results)
labels = {};
for ic = 1:numel(results)
    labels{end+1} = sprintf('%s Long', results(ic).name); %#ok<AGROW>
    labels{end+1} = sprintf('%s Short', results(ic).name); %#ok<AGROW>
end
end

function [ts, ys] = smooth_time_series(t, y)
t = t(:);
y = y(:);
if numel(t) < 4
    ts = t;
    ys = y;
    return;
end
ts = linspace(min(t), max(t), max(200, numel(t)*12)).';
ys = interp1(t, y, ts, 'pchip');
end

function fig = make_figure(sz)
if nargin < 1
    sz = [1000 700];
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, sz(1), sz(2)]);
set(fig, 'DefaultAxesFontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman');
set(fig, 'DefaultAxesFontSize', 16);
set(fig, 'DefaultTextFontSize', 16);
end

function save_fig(fig, pathName)
set(findall(fig, 'Type', 'axes'), 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'LineWidth', 1.1);
set(findall(fig, 'Type', 'text'), 'FontName', 'Times New Roman', 'FontWeight', 'bold');
exportgraphics(fig, pathName, 'Resolution', 220);
close(fig);
end

function write_metrics_csv(outDir, results)
pathName = fullfile(outDir, 'amplitude_comparison_metrics.csv');
fid = fopen(pathName, 'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'case,epsilon,epsilon_source,n_frames,late_start_frame,late_end_frame,dominant_frequency,mean_long_amp_f0,mean_short_amp_f0,long_short_amp_ratio,long_amp_over_epsilon,short_amp_over_epsilon,mean_long_amp_2f0,mean_short_amp_2f0,long_2f0_over_f0,short_2f0_over_f0,mean_long_rms,mean_short_rms\n');
for i = 1:numel(results)
    fprintf(fid, '%s,%.12g,%s,%d,%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n', ...
        results(i).name, results(i).epsilon, results(i).epsilonSource, results(i).nFrames, ...
        results(i).lateIdx(1), results(i).lateIdx(end), results(i).dominantFrequency, ...
        results(i).meanLongAmp, results(i).meanShortAmp, results(i).longShortAmpRatio, ...
        results(i).meanLongAmp/results(i).epsilon, results(i).meanShortAmp/results(i).epsilon, ...
        results(i).meanLongAmp2, results(i).meanShortAmp2, ...
        results(i).meanLongH2H1, results(i).meanShortH2H1, ...
        results(i).meanLongRms, results(i).meanShortRms);
end
end

function write_summary(outDir, results, f0, dtSave, fs)
pathName = fullfile(outDir, 'amplitude_comparison_summary.md');
fid = fopen(pathName, 'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '# E-2 / E-3 / E-4 amplitude-comparison summary\n\n');
fprintf(fid, 'Save interval: `%.12g`; sampling frequency: `%.12g`.\n\n', dtSave, fs);
fprintf(fid, 'The harmonic projection frequency is `f0 = %.12g`. E-2 and E-3 do not contain a local Config.cfg, so their epsilons are inferred from the folder names.\n\n', f0);
fprintf(fid, '| Case | epsilon | frames | late window | dominant FFT frequency | mean long | mean short | long/short | long/epsilon | short/epsilon | long 2f0/f0 | short 2f0/f0 |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i = 1:numel(results)
    fprintf(fid, '| %s | %.3g | %d | %d-%d | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g |\n', ...
        results(i).name, results(i).epsilon, results(i).nFrames, results(i).lateIdx(1), results(i).lateIdx(end), ...
        results(i).dominantFrequency, results(i).meanLongAmp, results(i).meanShortAmp, results(i).longShortAmpRatio, ...
        results(i).meanLongAmp/results(i).epsilon, results(i).meanShortAmp/results(i).epsilon, ...
        results(i).meanLongH2H1, results(i).meanShortH2H1);
end

epsVals = [results.epsilon];
[~, idxSmall] = min(epsVals);
fprintf(fid, '\n## Interpretation\n\n');
fprintf(fid, '- If the response is linear, `|P_f0|` should scale with `epsilon`, and `|P_f0|/epsilon` should remain nearly constant.\n');
fprintf(fid, '- The smallest-amplitude case `%s` is used as the linear reference for the nonlinearity-index plot.\n', results(idxSmall).name);
fprintf(fid, '- Growth of `|P_2f0|/|P_f0|` indicates waveform distortion and nonlinear generation of higher harmonics.\n');
fprintf(fid, '- Long/short-axis differences quantify how the elliptic cone geometry modulates the same incoming fast acoustic disturbance.\n\n');

for i = 1:numel(results)
    fprintf(fid, '## %s\n\n', results(i).name);
    fprintf(fid, 'Mean long-axis `|P_f0| = %.6g`, mean short-axis `|P_f0| = %.6g`, so the long/short ratio is `%.4g`.\n', ...
        results(i).meanLongAmp, results(i).meanShortAmp, results(i).longShortAmpRatio);
    fprintf(fid, 'The normalized responses are long `%.6g` and short `%.6g`.\n', ...
        results(i).meanLongAmp/results(i).epsilon, results(i).meanShortAmp/results(i).epsilon);
    fprintf(fid, 'The second-harmonic ratios are long `%.6g` and short `%.6g`.\n\n', ...
        results(i).meanLongH2H1, results(i).meanShortH2H1);
end
end
