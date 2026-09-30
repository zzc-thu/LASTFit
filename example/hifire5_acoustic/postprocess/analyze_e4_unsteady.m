function analyze_e4_unsteady()
% Analyze E-4 unsteady fast-wave case using RESU and Pert VTK outputs.

scriptDir = fileparts(mfilename('fullpath'));
caseReleaseDir = fileparts(scriptDir);
rawRoot = getenv('LASTFIT_HIFIRE_RAW_DIR');
if isempty(rawRoot)
    error(['Set LASTFIT_HIFIRE_RAW_DIR to the directory containing ', ...
        'E-2, E-3, E-4, and RESU (for example, Un_Hifire_fine).']);
end
rootDir = fullfile(rawRoot, 'E-4');
resuDir = fullfile(rootDir, 'RESU');
pertDir = fullfile(rootDir, 'Pert');
outDir = fullfile(caseReleaseDir, 'output', 'manuscript_fine', ...
    'analysis_unsteady');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

cfg = read_config(fullfile(rootDir, 'Config.cfg'));
pertFiles = list_pvts(pertDir, 'Pert_Result_(\d+)\.pvts');
resuFiles = list_pvts(resuDir, 'Result_(\d+)\.pvts');

fprintf('Pert frames: %d, RESU frames: %d\n', numel(pertFiles.idx), numel(resuFiles.idx));
fprintf('First/last Pert index: %d / %d\n', pertFiles.idx(1), pertFiles.idx(end));

firstInfo = parse_pvts(fullfile(pertDir, pertFiles.names{1}));
fields = firstInfo.fields;
disp('Pert fields:'); disp(strjoin(fields, ', '));
fieldP = choose_field(fields, {'P','P_pert','p','p_pert'});
fieldU = choose_field(fields, {'U','U_pert','u','u_pert'});
fieldV = choose_field(fields, {'V','V_pert','v','v_pert'});
fieldW = choose_field(fields, {'W','W_pert','w','w_pert'});
fieldT = choose_field(fields, {'T','T_pert','t','t_pert'});

dtOut = cfg.StepsWriteData * cfg.FixedDeltaTime;
if ~isfinite(dtOut) || dtOut <= 0
    dtOut = 1;
end
t = (0:numel(pertFiles.idx)-1)' * dtOut;

% Axis indices identified from the steady wall geometry in the prior analysis.
kLabels = {'Long axis +','Long axis -','Short axis +','Short axis -'};
kList = [41 121 1 81];
iStations = [40 81 121];
wallJ = 1;

nt = numel(pertFiles.idx);
ni = firstInfo.whole(2) - firstInfo.whole(1) + 1;
nk = firstInfo.whole(6) - firstInfo.whole(5) + 1;
Pline = nan(nt, ni, numel(kList));
UmagLine = nan(nt, ni, numel(kList));
Tline = nan(nt, ni, numel(kList));
Xline = nan(ni, numel(kList));

lineFields = {fieldP, fieldU, fieldV, fieldW, fieldT};
for it = 1:nt
    pvtsFile = fullfile(pertDir, pertFiles.names{it});
    for kk = 1:numel(kList)
        line = read_wall_line(pvtsFile, pertDir, kList(kk), wallJ, lineFields);
        Xline(:,kk) = line.X(:);
        Pline(it,:,kk) = line.(fieldP)(:);
        Tline(it,:,kk) = line.(fieldT)(:);
        UmagLine(it,:,kk) = sqrt(line.(fieldU)(:).^2 + line.(fieldV)(:).^2 + line.(fieldW)(:).^2);
    end
    if mod(it, max(1,round(nt/10))) == 0 || it == nt
        fprintf('  line time series %d/%d\n', it, nt);
    end
end

% RESU total-pressure reference at the last common frame for relative amplitudes.
[hasRef, refP] = read_resu_reference(resuFiles, pertFiles.idx(end), resuDir, kList, wallJ, fieldP);

Prms = squeeze(sqrt(mean(Pline.^2, 1, 'omitnan')));
Pmax = squeeze(max(abs(Pline), [], 1, 'omitnan'));
Urms = squeeze(sqrt(mean(UmagLine.^2, 1, 'omitnan')));
Trms = squeeze(sqrt(mean(Tline.^2, 1, 'omitnan')));

% Temporal spectral content at representative streamwise station.
iFFT = iStations(2);
sig = squeeze(Pline(:,iFFT,1));
[domFreq, domAmp, phaseLine] = dominant_frequency_and_phase(Pline, t, iFFT);

% Sampled full-wall RMS map of perturbation pressure.
nMap = min(16, nt);
mapIds = unique(round(linspace(1, nt, nMap)));
PwallSumSq = zeros(ni, nk);
PwallCount = 0;
Xwall = nan(ni, nk);
for jj = 1:numel(mapIds)
    it = mapIds(jj);
    pvtsFile = fullfile(pertDir, pertFiles.names{it});
    [Xtmp, Ptmp] = read_wall_map(pvtsFile, pertDir, fieldP, wallJ);
    Xwall = Xtmp;
    PwallSumSq = PwallSumSq + Ptmp.^2;
    PwallCount = PwallCount + 1;
    fprintf('  wall RMS sample %d/%d\n', jj, numel(mapIds));
end
PwallRms = sqrt(PwallSumSq / max(1,PwallCount));

write_summary(outDir, cfg, pertFiles, resuFiles, kLabels, kList, Xline, Prms, Pmax, Urms, Trms, ...
    hasRef, refP, iStations, Pline, t, domFreq, domAmp, phaseLine, PwallRms);
write_csv(outDir, Xline, kLabels, Prms, Pmax, Urms, Trms, hasRef, refP);

plot_axis_rms(outDir, Xline, kLabels, Prms, 'P perturbation RMS', 'P''_{rms}', 'axis_pressure_rms.png');
plot_axis_rms(outDir, Xline, kLabels, Pmax, 'P perturbation peak amplitude', 'max |P''|', 'axis_pressure_peak.png');
plot_axis_rms(outDir, Xline, kLabels, Urms, 'Velocity perturbation RMS', '|u''|_{rms}', 'axis_velocity_rms.png');
plot_time_histories(outDir, t, Pline, Xline, kLabels, iStations);
plot_phase(outDir, Xline, kLabels, phaseLine, domFreq);
plot_wall_rms_map(outDir, Xwall, PwallRms, kList, kLabels);

fprintf('Done. Outputs in %s\n', outDir);
end

function cfg = read_config(fn)
txt = fileread(fn);
cfg.Pert_Type = getnum(txt, 'Pert_Type');
cfg.k_infty = getnum(txt, 'k_infty');
cfg.epsilon = getnum(txt, 'epsilon');
cfg.Mach_Ref = getnum(txt, 'Mach_Ref');
cfg.Re_Ref = getnum(txt, 'Re_Ref');
cfg.TWall = getnum(txt, 'TWall');
cfg.FixedDeltaTime = getnum(txt, 'FixedDeltaTime');
cfg.StepsWriteData = getnum(txt, 'StepsWriteData');
cfg.Nx = getnum(txt, 'Nx');
cfg.Ny = getnum(txt, 'Ny');
cfg.Nz = getnum(txt, 'Nz');
end

function val = getnum(txt, key)
expr = ['(^|\n)\s*' key '\s*=\s*([+-]?\d*\.?\d+(?:[dDeE][+-]?\d+)?)'];
tok = regexp(txt, expr, 'tokens', 'once');
if isempty(tok)
    val = NaN;
else
    val = str2double(regexprep(tok{2}, '[dD]', 'e'));
end
end

function files = list_pvts(folder, expr)
d = dir(fullfile(folder, '*.pvts'));
idx = [];
names = {};
for i = 1:numel(d)
    tok = regexp(d(i).name, expr, 'tokens', 'once');
    if ~isempty(tok)
        idx(end+1,1) = str2double(tok{1}); %#ok<AGROW>
        names{end+1,1} = d(i).name; %#ok<AGROW>
    end
end
[idx, order] = sort(idx);
files.idx = idx;
files.names = names(order);
end

function name = choose_field(fields, candidates)
name = '';
for c = 1:numel(candidates)
    hit = strcmp(fields, candidates{c});
    if any(hit)
        name = fields{find(hit,1)};
        return;
    end
end
error('Cannot find field among candidates: %s', strjoin(candidates, ', '));
end

function info = parse_pvts(pvtsFile)
txt = fileread(pvtsFile);
tok = regexp(txt, 'WholeExtent="\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)', 'tokens', 'once');
info.whole = str2double(tok);
ptoks = regexp(txt, '<Piece\s+Extent=\s*"\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+"\s+Source="([^"]+)"', 'tokens');
pieces = repmat(struct('extent', [], 'source', ''), 1, numel(ptoks));
for n = 1:numel(ptoks)
    pieces(n).extent = str2double(ptoks{n}(1:6));
    pieces(n).source = ptoks{n}{7};
end
info.pieces = pieces;
ftoks = regexp(txt, '<PDataArray[^>]*Name="([^"]+)"', 'tokens');
info.fields = cellfun(@(x)x{1}, ftoks, 'UniformOutput', false);
end

function line = read_wall_line(pvtsFile, folder, kIndex, wallJ, fields)
info = parse_pvts(pvtsFile);
ni = info.whole(2) - info.whole(1) + 1;
globalK = kIndex - 1 + info.whole(5);
globalJ = wallJ - 1 + info.whole(3);
line.X = nan(ni,1);
for f = 1:numel(fields)
    line.(fields{f}) = nan(ni,1);
end
for p = 1:numel(info.pieces)
    e = info.pieces(p).extent;
    if globalK < e(5) || globalK > e(6) || globalJ < e(3) || globalJ > e(4)
        continue;
    end
    data = read_vts_piece(fullfile(folder, info.pieces(p).source), fields);
    iiAll = (e(1):e(2)) - info.whole(1) + 1;
    validI = iiAll >= 1 & iiAll <= ni;
    ii = iiAll(validI);
    jj = globalJ - e(3) + 1;
    kk = globalK - e(5) + 1;
    line.X(ii) = double(squeeze(data.X(validI,jj,kk)));
    for f = 1:numel(fields)
        vals = squeeze(data.(fields{f})(validI,jj,kk));
        line.(fields{f})(ii) = double(vals);
    end
end
end

function [Xwall, Pwall] = read_wall_map(pvtsFile, folder, fieldP, wallJ)
info = parse_pvts(pvtsFile);
ni = info.whole(2) - info.whole(1) + 1;
nk = info.whole(6) - info.whole(5) + 1;
globalJ = wallJ - 1 + info.whole(3);
Xwall = nan(ni,nk);
Pwall = nan(ni,nk);
for p = 1:numel(info.pieces)
    e = info.pieces(p).extent;
    if globalJ < e(3) || globalJ > e(4)
        continue;
    end
    data = read_vts_piece(fullfile(folder, info.pieces(p).source), {fieldP});
    iiAll = (e(1):e(2)) - info.whole(1) + 1;
    kkAll = (e(5):e(6)) - info.whole(5) + 1;
    validI = iiAll >= 1 & iiAll <= ni;
    validK = kkAll >= 1 & kkAll <= nk;
    ii = iiAll(validI);
    kk = kkAll(validK);
    jj = globalJ - e(3) + 1;
    Xwall(ii,kk) = double(squeeze(data.X(validI,jj,validK)));
    Pwall(ii,kk) = double(squeeze(data.(fieldP)(validI,jj,validK)));
end
end

function data = read_vts_piece(fn, fields)
fid = fopen(fn, 'r', 'ieee-le');
if fid < 0
    error('Cannot open %s', fn);
end
bytes = fread(fid, inf, 'uint8=>uint8');
fclose(fid);
pat = uint8('<AppendedData encoding="raw">_');
idx = strfind(bytes', pat);
if isempty(idx)
    error('No appended data in %s', fn);
end
base = idx(1) + numel(pat);
header = char(bytes(1:idx(1)-1))';
extTok = regexp(header, 'Piece\s+Extent="\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)', 'tokens', 'once');
extent = str2double(extTok);
ni = extent(2) - extent(1) + 1;
nj = extent(4) - extent(3) + 1;
nk = extent(6) - extent(5) + 1;
npts = ni * nj * nk;
ptTok = regexp(header, '<Points>\s*<DataArray[^>]*offset="\s*(\d+)', 'tokens', 'once');
ptOffset = str2double(ptTok{1});
ptStart = base + ptOffset + 8;
pts = typecast(bytes(ptStart:ptStart + 3*npts*4 - 1), 'single');
pts = reshape(pts, 3, npts);
data.X = reshape(pts(1,:), ni, nj, nk);
for i = 1:numel(fields)
    nm = fields{i};
    expr = ['<DataArray[^>]*Name="' nm '"[^>]*offset="\s*(\d+)'];
    tok = regexp(header, expr, 'tokens', 'once');
    if isempty(tok)
        error('Field %s missing in %s', nm, fn);
    end
    off = str2double(tok{1});
    start = base + off + 8;
    vals = typecast(bytes(start:start + npts*4 - 1), 'single');
    data.(nm) = reshape(vals, ni, nj, nk);
end
end

function [hasRef, refP] = read_resu_reference(resuFiles, targetIdx, resuDir, kList, wallJ, fieldP)
hasRef = false;
refP = [];
if isempty(resuFiles.idx)
    return;
end
[~, loc] = min(abs(resuFiles.idx - targetIdx));
try
    pvtsFile = fullfile(resuDir, resuFiles.names{loc});
    info = parse_pvts(pvtsFile);
    ni = info.whole(2) - info.whole(1) + 1;
    refP = nan(ni, numel(kList));
    for kk = 1:numel(kList)
        line = read_wall_line(pvtsFile, resuDir, kList(kk), wallJ, {fieldP});
        refP(:,kk) = line.(fieldP);
    end
    hasRef = true;
catch ME
    warning('RESU reference read failed: %s', ME.message);
end
end

function [domFreq, domAmp, phaseLine] = dominant_frequency_and_phase(Pline, t, iFFT)
nt = size(Pline,1);
dt = median(diff(t));
sig = squeeze(Pline(:,iFFT,1));
sig = sig - mean(sig,'omitnan');
Y = fft(sig);
freq = (0:nt-1)'/(nt*dt);
half = 2:floor(nt/2);
[domAmp, loc] = max(abs(Y(half))/nt);
domIdx = half(loc);
domFreq = freq(domIdx);
phaseLine = nan(size(Pline,2), size(Pline,3));
for kk = 1:size(Pline,3)
    for i = 1:size(Pline,2)
        s = squeeze(Pline(:,i,kk));
        s = s - mean(s,'omitnan');
        F = fft(s);
        phaseLine(i,kk) = angle(F(domIdx));
    end
end
end

function write_summary(outDir, cfg, pertFiles, resuFiles, labels, kList, Xline, Prms, Pmax, Urms, Trms, hasRef, refP, iStations, Pline, t, domFreq, domAmp, phaseLine, PwallRms)
fid = fopen(fullfile(outDir, 'e4_unsteady_summary.txt'), 'w');
fprintf(fid, 'E-4 unsteady fast-wave analysis\n\n');
fprintf(fid, 'Config: Pert_Type=%g, k_infty=%g, epsilon=%g, Mach=%g, Re=%g, TWall=%g\n', ...
    cfg.Pert_Type, cfg.k_infty, cfg.epsilon, cfg.Mach_Ref, cfg.Re_Ref, cfg.TWall);
fprintf(fid, 'Time output: FixedDeltaTime=%g, StepsWriteData=%g, assumed dt_between_frames=%g\n', ...
    cfg.FixedDeltaTime, cfg.StepsWriteData, cfg.FixedDeltaTime*cfg.StepsWriteData);
fprintf(fid, 'Frames: Pert=%d [%d..%d], RESU=%d [%d..%d]\n\n', ...
    numel(pertFiles.idx), pertFiles.idx(1), pertFiles.idx(end), numel(resuFiles.idx), resuFiles.idx(1), resuFiles.idx(end));
fprintf(fid, 'Dominant temporal frequency at i=%d, %s: f=%.8g, FFT amplitude=%.8g\n\n', ...
    iStations(2), labels{1}, domFreq, domAmp);

validI = 10:size(Prms,1);
fprintf(fid, 'Wall-meridian perturbation amplitudes, excluding first 9 i-points:\n');
for kk = 1:numel(labels)
    fprintf(fid, '  %-12s k=%3d mean(P_rms)=%.8g max(P_rms)=%.8g max(|P|)=%.8g mean(|u|_rms)=%.8g mean(T_rms)=%.8g\n', ...
        labels{kk}, kList(kk), mean(Prms(validI,kk),'omitnan'), max(Prms(validI,kk)), ...
        max(Pmax(validI,kk)), mean(Urms(validI,kk),'omitnan'), mean(Trms(validI,kk),'omitnan'));
end
longAvg = mean(Prms(validI,1:2), 2, 'omitnan');
shortAvg = mean(Prms(validI,3:4), 2, 'omitnan');
fprintf(fid, '\nLong/short pressure RMS anisotropy:\n');
fprintf(fid, '  mean(P_rms_long/P_rms_short)=%.8g\n', mean(longAvg./shortAvg,'omitnan'));
fprintf(fid, '  mean(P_rms_long-P_rms_short)=%.8g\n', mean(longAvg-shortAvg,'omitnan'));
if hasRef
    refMean = mean(abs(refP(validI,:)),1,'omitnan');
    rel = mean(Prms(validI,:),1,'omitnan') ./ refMean;
    fprintf(fid, '\nRelative to RESU total wall pressure at reference frame:\n');
    for kk = 1:numel(labels)
        fprintf(fid, '  %-12s mean(P_rms)/mean(|P_total|)=%.8g\n', labels{kk}, rel(kk));
    end
end
fprintf(fid, '\nSample station time-series peak-to-peak P perturbation:\n');
for ii = iStations
    for kk = [1 3]
        s = squeeze(Pline(:,ii,kk));
        fprintf(fid, '  i=%3d %-12s X=%.8g peak_to_peak=%.8g rms=%.8g\n', ...
            ii, labels{kk}, Xline(ii,kk), max(s)-min(s), sqrt(mean(s.^2,'omitnan')));
    end
end
fprintf(fid, '\nWall RMS map: max(P_rms)=%.8g, mean(P_rms)=%.8g\n', max(PwallRms(:)), mean(PwallRms(:),'omitnan'));
fclose(fid);
end

function write_csv(outDir, Xline, labels, Prms, Pmax, Urms, Trms, hasRef, refP)
fid = fopen(fullfile(outDir, 'axis_perturbation_profiles.csv'), 'w');
fprintf(fid, 'i,x,label,P_rms,P_peak,Umag_rms,T_rms');
if hasRef
    fprintf(fid, ',P_total_ref,relative_P_rms');
end
fprintf(fid, '\n');
for kk = 1:numel(labels)
    for i = 1:size(Prms,1)
        fprintf(fid, '%d,%.10g,%s,%.10g,%.10g,%.10g,%.10g', ...
            i, Xline(i,kk), labels{kk}, Prms(i,kk), Pmax(i,kk), Urms(i,kk), Trms(i,kk));
        if hasRef
            fprintf(fid, ',%.10g,%.10g', refP(i,kk), Prms(i,kk)/abs(refP(i,kk)));
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);
end

function plot_axis_rms(outDir, Xline, labels, A, ttl, ylab, fname)
fig = figure('Visible','off','Position',[100,100,1100,720]);
styles = {'-','-','--','--'};
for kk = 1:numel(labels)
    plot(Xline(:,kk), A(:,kk), styles{kk}, 'LineWidth', 2.2); hold on;
end
xlabel('X'); ylabel(ylab); title(ttl);
legend(labels, 'Location','best');
grid on; box on; style_pub(gca);
print(fig, fullfile(outDir, fname), '-dpng', '-r220');
close(fig);
end

function plot_time_histories(outDir, t, Pline, Xline, labels, iStations)
fig = figure('Visible','off','Position',[100,100,1800,1180], 'Color','w');
legendHandles = [];
legendLabels = {labels{1}, labels{3}};
legendAx = [];
axPos = [
    0.135 0.650 0.835 0.205
    0.135 0.395 0.835 0.205
    0.135 0.140 0.835 0.205
];
for s = 1:numel(iStations)
    ii = iStations(s);
    ax = axes('Parent', fig, 'Position', axPos(s,:)); %#ok<LAXES>
    if s == 1
        legendAx = ax;
    end
    yLong = squeeze(Pline(:,ii,1));
    yShort = squeeze(Pline(:,ii,3));
    [ts, yLongS] = hybrid_periodic_time_series(t, yLong);
    [~, yShortS] = hybrid_periodic_time_series(t, yShort);
    yScale = 10^ceil(-log10(max(abs([yLong(:); yShort(:)]), [], 'omitnan')));
    if ~isfinite(yScale) || yScale <= 0
        yScale = 1;
    end
    yExp = -round(log10(yScale));
    hLong = plot(ts, yLongS*yScale, '-', 'LineWidth', 3.4); hold on;
    hShort = plot(ts, yShortS*yScale, '--', 'LineWidth', 3.4);
    if s == 1
        legendHandles = [hLong, hShort]; %#ok<AGROW>
    else
        set([hLong, hShort], 'HandleVisibility','off');
    end
    plot(t, yLong*yScale, 'o', 'MarkerSize', 3.2, 'LineWidth', 0.8, 'HandleVisibility','off');
    plot(t, yShort*yScale, 's', 'MarkerSize', 3.2, 'LineWidth', 0.8, 'HandleVisibility','off');
    ylabel(sprintf('P'' (10^{%d})', yExp), 'Interpreter','tex');
    text(0.025, 0.78, sprintf('\\iti\\rm = %d', ii), 'Units','normalized', ...
        'Interpreter','tex', 'FontName','Times New Roman', ...
        'FontSize',22, 'FontWeight','normal', 'BackgroundColor','w', 'Margin',2);
    grid on; box on;
    set(ax, 'FontName','Times New Roman', 'FontSize',20, 'FontWeight','normal', 'LineWidth',1.1);
    set(get(ax,'YLabel'), 'FontName','Times New Roman', 'FontSize',24, 'FontWeight','normal');
    if s ~= numel(iStations)
        set(ax, 'XTickLabel', []);
    end
end
lgd = legend(legendAx, legendHandles, legendLabels, 'Orientation','horizontal', 'Location','northoutside');
set(lgd, 'FontName','Times New Roman','FontSize',24,'FontWeight','normal', ...
    'Box','on', 'Units','normalized', 'Position',[0.37 0.860 0.26 0.048]);
annotation(fig, 'textbox', [0.06 0.925 0.88 0.045], ...
    'String','Wall pressure perturbation time history', ...
    'FontName','Times New Roman', 'FontSize',30, 'FontWeight','normal', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'LineStyle','none');
annotation(fig, 'textbox', [0.43 0.030 0.14 0.045], ...
    'String','time', ...
    'FontName','Times New Roman', 'FontSize',28, 'FontWeight','normal', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'LineStyle','none');
print(fig, fullfile(outDir, 'pressure_perturbation_time_histories.png'), '-dpng', '-r220');
close(fig);

plot_time_histories_raw(outDir, t, Pline, labels, iStations);
end

function plot_time_histories_raw(outDir, t, Pline, labels, iStations)
fig = figure('Visible','off','Position',[100,100,1800,1180], 'Color','w');
legendHandles = [];
legendLabels = {labels{1}, labels{3}};
legendAx = [];
axPos = [
    0.135 0.650 0.835 0.205
    0.135 0.395 0.835 0.205
    0.135 0.140 0.835 0.205
];
for s = 1:numel(iStations)
    ii = iStations(s);
    ax = axes('Parent', fig, 'Position', axPos(s,:)); %#ok<LAXES>
    if s == 1
        legendAx = ax;
    end
    yLong = squeeze(Pline(:,ii,1));
    yShort = squeeze(Pline(:,ii,3));
    yScale = 10^ceil(-log10(max(abs([yLong(:); yShort(:)]), [], 'omitnan')));
    if ~isfinite(yScale) || yScale <= 0
        yScale = 1;
    end
    yExp = -round(log10(yScale));
    hLong = plot(t, yLong*yScale, 'o', 'MarkerSize', 5.2, ...
        'LineWidth', 1.2, 'MarkerFaceColor','none'); hold on;
    hShort = plot(t, yShort*yScale, 's', 'MarkerSize', 5.2, ...
        'LineWidth', 1.2, 'MarkerFaceColor','none');
    if s == 1
        legendHandles = [hLong, hShort]; %#ok<AGROW>
    else
        set([hLong, hShort], 'HandleVisibility','off');
    end
    ylabel(sprintf('P'' (10^{%d})', yExp), 'Interpreter','tex');
    text(0.025, 0.78, sprintf('\\iti\\rm = %d', ii), 'Units','normalized', ...
        'Interpreter','tex', 'FontName','Times New Roman', ...
        'FontSize',22, 'FontWeight','normal', 'BackgroundColor','w', 'Margin',2);
    grid on; box on;
    set(ax, 'FontName','Times New Roman', 'FontSize',20, ...
        'FontWeight','normal', 'LineWidth',1.1);
    set(get(ax,'YLabel'), 'FontName','Times New Roman', ...
        'FontSize',24, 'FontWeight','normal');
    if s ~= numel(iStations)
        set(ax, 'XTickLabel', []);
    end
end
lgd = legend(legendAx, legendHandles, legendLabels, ...
    'Orientation','horizontal', 'Location','northoutside');
set(lgd, 'FontName','Times New Roman','FontSize',24,'FontWeight','normal', ...
    'Box','on', 'Units','normalized', 'Position',[0.37 0.860 0.26 0.048]);
annotation(fig, 'textbox', [0.06 0.925 0.88 0.045], ...
    'String','Wall pressure perturbation time history (raw data)', ...
    'FontName','Times New Roman', 'FontSize',30, 'FontWeight','normal', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'LineStyle','none');
annotation(fig, 'textbox', [0.43 0.030 0.14 0.045], ...
    'String','time', ...
    'FontName','Times New Roman', 'FontSize',28, 'FontWeight','normal', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'LineStyle','none');
print(fig, fullfile(outDir, 'pressure_perturbation_time_histories_raw.png'), '-dpng', '-r220');
close(fig);
end

function [ts, ys] = smooth_time_series(t, y)
ok = isfinite(t(:)) & isfinite(y(:));
t0 = t(ok);
y0 = y(ok);
if numel(t0) < 4
    ts = t0;
    ys = y0;
    return;
end
ts = linspace(min(t0), max(t0), max(1600, numel(t0)*60))';
ys0 = interp1(t0, y0, ts, 'makima');
window = max(9, round(numel(ts)/180));
if mod(window, 2) == 0
    window = window + 1;
end
ys = smoothdata(ys0, 'sgolay', window);
end

function [ts, ys] = hybrid_periodic_time_series(t, y)
ok = isfinite(t(:)) & isfinite(y(:));
t0 = t(ok);
y0 = y(ok);
if numel(t0) < 4
    ts = t0;
    ys = y0;
    return;
end

ts = linspace(min(t0), max(t0), max(1800, numel(t0)*60))';
ysRaw = interp1(t0, y0, ts, 'pchip');
ysRaw = min(max(ysRaw, min(y0)), max(y0));

ampRef = max(abs(y0));
active = abs(y0) > max(0.08*ampRef, eps);
firstActive = find(active, 1, 'first');
f0 = 666.6666666667;
period = 1/f0;
if isempty(firstActive)
    ys = ysRaw;
    return;
end

% Keep the startup part raw, then fit the later nearly periodic response.
tFitStart = t0(firstActive) + 0.75*period;
fitIdx = t0 >= tFitStart;
if nnz(fitIdx) < 8
    ys = ysRaw;
    return;
end

omega = 2*pi*f0;
tf = t0(fitIdx);
yf = y0(fitIdx);
H = [ones(numel(tf),1), cos(omega*tf), sin(omega*tf), ...
     cos(2*omega*tf), sin(2*omega*tf)];
coef = H \ yf;
ysFit = [ones(numel(ts),1), cos(omega*ts), sin(omega*ts), ...
         cos(2*omega*ts), sin(2*omega*ts)] * coef;

% Match fitted amplitude to the raw data range in the fitted segment.
rawMin = min(yf);
rawMax = max(yf);
fitMaskDense = ts >= tFitStart;
fitMin = min(ysFit(fitMaskDense));
fitMax = max(ysFit(fitMaskDense));
if fitMax > fitMin && rawMax > rawMin
    ysFit = (ysFit - fitMin) / (fitMax - fitMin) * (rawMax - rawMin) + rawMin;
end
ysFit = min(max(ysFit, min(y0)), max(y0));

blendWidth = 0.35*period;
w = min(max((ts - tFitStart) / blendWidth, 0), 1);
w = w.^2 .* (3 - 2*w);
ys = (1-w).*ysRaw + w.*ysFit;
ys = min(max(ys, min(y0)), max(y0));
end

function plot_phase(outDir, Xline, labels, phaseLine, domFreq)
fig = figure('Visible','off','Position',[100,100,1100,720]);
for kk = [1 3]
    ph = unwrap(phaseLine(:,kk));
    plot(Xline(:,kk), ph, 'LineWidth', 2.2); hold on;
end
xlabel('X'); ylabel('phase (rad)');
title(sprintf('Dominant-frequency phase, f = %.4g', domFreq));
legend({labels{1}, labels{3}}, 'Location','best');
grid on; box on; style_pub(gca);
print(fig, fullfile(outDir, 'pressure_perturbation_phase.png'), '-dpng', '-r220');
close(fig);
end

function plot_wall_rms_map(outDir, Xwall, PwallRms, kList, labels)
xLine = mean(Xwall,2,'omitnan');
fig = figure('Visible','off','Position',[100,100,1300,820]);
imagesc(xLine, 1:size(PwallRms,2), PwallRms');
axis xy; hold on;
lineStyles = {'-k','-k','--k','--k'};
for kk = 1:numel(kList)
    plot(xLine, kList(kk)*ones(size(xLine)), lineStyles{kk}, 'LineWidth', 1.6);
    text(max(xLine)*0.985, kList(kk)+2, labels{kk}, 'HorizontalAlignment','right', ...
        'FontName','Times New Roman','FontSize',16,'FontWeight','bold','BackgroundColor','w');
end
xlabel('X'); ylabel('circumferential index k'); title('Wall pressure perturbation RMS');
cb = colorbar; ylabel(cb, 'P''_{rms}');
grid on; box on; style_pub(gca);
set(cb, 'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
print(fig, fullfile(outDir, 'wall_pressure_perturbation_rms_map.png'), '-dpng', '-r220');
close(fig);
end

function style_pub(ax)
set(ax, 'FontName','Times New Roman', 'FontSize',20, 'FontWeight','bold', 'LineWidth',1.3);
set(get(ax,'XLabel'), 'FontName','Times New Roman', 'FontSize',24, 'FontWeight','bold');
set(get(ax,'YLabel'), 'FontName','Times New Roman', 'FontSize',24, 'FontWeight','bold');
set(get(ax,'Title'), 'FontName','Times New Roman', 'FontSize',24, 'FontWeight','bold');
lgd = findobj(ancestor(ax,'figure'), 'Type','Legend');
for i = 1:numel(lgd)
    set(lgd(i), 'FontName','Times New Roman', 'FontSize',14, 'FontWeight','bold');
end
end
