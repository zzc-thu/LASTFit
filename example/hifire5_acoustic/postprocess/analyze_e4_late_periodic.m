function analyze_e4_late_periodic()
% Late-window periodic analysis for E-4 fast acoustic-wave perturbation.

scriptDir = fileparts(mfilename('fullpath'));
caseReleaseDir = fileparts(scriptDir);
rawRoot = getenv('LASTFIT_HIFIRE_RAW_DIR');
if isempty(rawRoot)
    error(['Set LASTFIT_HIFIRE_RAW_DIR to the directory containing ', ...
        'E-2, E-3, E-4, and RESU (for example, Un_Hifire_fine).']);
end
rootDir = fullfile(rawRoot, 'E-4');
pertDir = fullfile(rootDir, 'Pert');
outDir = fullfile(caseReleaseDir, 'output', 'manuscript_fine', ...
    'analysis_late_periodic');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

cfg = read_config(fullfile(rootDir, 'Config.cfg'));
files = list_pvts(pertDir, 'Pert_Result_(\d+)\.pvts');
nt = numel(files.idx);
dtOut = cfg.FixedDeltaTime * cfg.StepsWriteData;
t = (0:nt-1)' * dtOut;
lateIds = (floor(nt/2)+1):nt;
tLate = t(lateIds);

info = parse_pvts(fullfile(pertDir, files.names{1}));
ni = info.whole(2) - info.whole(1) + 1;
nk = info.whole(6) - info.whole(5) + 1;
fieldP = 'P_pert';
wallJ = 1;

kLabels = {'Long axis +','Long axis -','Short axis +','Short axis -'};
kList = [41 121 1 81];
iStations = [40 81 121];

fprintf('Late periodic analysis: total frames=%d, late frames=%d..%d (%d frames)\n', ...
    nt, lateIds(1), lateIds(end), numel(lateIds));

% Line time series on long/short meridians.
Pline = nan(nt, ni, numel(kList));
Xline = nan(ni, numel(kList));
for it = 1:nt
    pvtsFile = fullfile(pertDir, files.names{it});
    for kk = 1:numel(kList)
        line = read_wall_line(pvtsFile, pertDir, kList(kk), wallJ, {fieldP});
        Xline(:,kk) = line.X(:);
        Pline(it,:,kk) = line.(fieldP)(:);
    end
    if mod(it, 5) == 0 || it == nt
        fprintf('  meridian time series %d/%d\n', it, nt);
    end
end

% Dominant frequency from the full response, then harmonic projection in late window.
iFreq = iStations(2);
f0 = dominant_frequency(squeeze(Pline(:,iFreq,1)), t);
fprintf('Dominant full-window frequency f0 = %.8g\n', f0);

PmeanLate = squeeze(mean(Pline(lateIds,:,:), 1, 'omitnan'));
PrmsLate = squeeze(sqrt(mean(Pline(lateIds,:,:).^2, 1, 'omitnan')));
[PampLate, PphaseLate] = harmonic_amp_phase(Pline(lateIds,:,:), tLate, f0);

% Wall harmonic maps over the late window.
Psum = zeros(ni,nk);
PsumSq = zeros(ni,nk);
Pcomplex = zeros(ni,nk);
Xwall = nan(ni,nk);
omega = 2*pi*f0;
for jj = 1:numel(lateIds)
    it = lateIds(jj);
    pvtsFile = fullfile(pertDir, files.names{it});
    [Xtmp, Ptmp] = read_wall_map(pvtsFile, pertDir, fieldP, wallJ);
    Xwall = Xtmp;
    Psum = Psum + Ptmp;
    PsumSq = PsumSq + Ptmp.^2;
    Pcomplex = Pcomplex + Ptmp .* exp(-1i * omega * t(it));
    fprintf('  late wall map %d/%d\n', jj, numel(lateIds));
end
Nlate = numel(lateIds);
PwallMeanLate = Psum / Nlate;
PwallRmsLate = sqrt(PsumSq / Nlate);
PwallAmpLate = 2 * abs(Pcomplex) / Nlate;
PwallPhaseLate = angle(Pcomplex);

% Window-to-window RMS convergence over late half.
conv = late_window_convergence(Pline, lateIds, iStations, kLabels);

write_late_summary(outDir, cfg, files, lateIds, tLate, f0, kLabels, kList, ...
    Xline, PrmsLate, PampLate, PphaseLate, iStations, conv, PwallRmsLate, PwallAmpLate);
write_late_csv(outDir, Xline, kLabels, PmeanLate, PrmsLate, PampLate, PphaseLate);

plot_axis_quantity(outDir, Xline, kLabels, PrmsLate, 'Late-window wall pressure perturbation RMS', 'P''_{rms}', 'late_axis_pressure_rms.png');
plot_axis_quantity(outDir, Xline, kLabels, PampLate, 'Late-window harmonic amplitude of wall pressure perturbation', '|P''_{f0}|', 'late_axis_pressure_harmonic_amplitude.png');
plot_axis_phase(outDir, Xline, kLabels, PphaseLate, f0);
plot_wall_map(outDir, Xwall, PwallRmsLate, kList, kLabels, 'Late-window wall pressure perturbation RMS', 'P''_{rms}', 'late_wall_pressure_rms_map.png');
plot_wall_map(outDir, Xwall, PwallAmpLate, kList, kLabels, 'Late-window harmonic amplitude of wall pressure perturbation', '|P''_{f0}|', 'late_wall_pressure_harmonic_amplitude_map.png');
plot_wall_phase(outDir, Xwall, PwallPhaseLate, kList, kLabels, f0);
plot_late_time_fit(outDir, t, Pline, Xline, kLabels, iStations, lateIds, f0);
plot_convergence(outDir, conv, iStations);

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

function f0 = dominant_frequency(sig, t)
sig = sig(:) - mean(sig,'omitnan');
nt = numel(sig);
dt = median(diff(t));
Y = fft(sig);
freq = (0:nt-1)'/(nt*dt);
half = 2:floor(nt/2);
[~, loc] = max(abs(Y(half)));
f0 = freq(half(loc));
end

function [amp, phase] = harmonic_amp_phase(A, t, f0)
omega = 2*pi*f0;
N = numel(t);
shape = size(A);
amp = nan(shape(2), shape(3));
phase = nan(shape(2), shape(3));
for kk = 1:shape(3)
    for i = 1:shape(2)
        sig = squeeze(A(:,i,kk));
        sig = sig - mean(sig,'omitnan');
        c = sum(sig .* exp(-1i*omega*t(:)), 'omitnan');
        amp(i,kk) = 2 * abs(c) / N;
        phase(i,kk) = angle(c);
    end
end
end

function conv = late_window_convergence(Pline, lateIds, iStations, labels)
nWin = 3;
edges = round(linspace(1, numel(lateIds)+1, nWin+1));
conv.rms = nan(nWin, numel(iStations), 2);
conv.labels = {labels{1}, labels{3}};
for w = 1:nWin
    ids = lateIds(edges(w):edges(w+1)-1);
    for s = 1:numel(iStations)
        ii = iStations(s);
        conv.rms(w,s,1) = sqrt(mean(squeeze(Pline(ids,ii,1)).^2, 'omitnan'));
        conv.rms(w,s,2) = sqrt(mean(squeeze(Pline(ids,ii,3)).^2, 'omitnan'));
    end
end
end

function write_late_summary(outDir, cfg, files, lateIds, tLate, f0, labels, kList, Xline, Prms, Pamp, Pphase, iStations, conv, PwallRms, PwallAmp)
fid = fopen(fullfile(outDir, 'late_periodic_summary.txt'), 'w');
fprintf(fid, 'E-4 late-window periodic acoustic perturbation analysis\n\n');
fprintf(fid, 'Config: Pert_Type=%g fast wave, k_infty=%g, epsilon=%g, Mach=%g, Re=%g, TWall=%g\n', ...
    cfg.Pert_Type, cfg.k_infty, cfg.epsilon, cfg.Mach_Ref, cfg.Re_Ref, cfg.TWall);
fprintf(fid, 'Output dt=%g, total frames=%d, late frames=%d..%d, late time=[%.8g, %.8g]\n', ...
    cfg.FixedDeltaTime*cfg.StepsWriteData, numel(files.idx), lateIds(1), lateIds(end), tLate(1), tLate(end));
fprintf(fid, 'Dominant full-window frequency used for late harmonic projection: f0=%.8g\n\n', f0);
validI = 10:size(Prms,1);
for kk = 1:numel(labels)
    fprintf(fid, '%-12s k=%3d mean(P_rms)=%.8g max(P_rms)=%.8g mean(|P_f0|)=%.8g max(|P_f0|)=%.8g\n', ...
        labels{kk}, kList(kk), mean(Prms(validI,kk),'omitnan'), max(Prms(validI,kk)), ...
        mean(Pamp(validI,kk),'omitnan'), max(Pamp(validI,kk)));
end
longAmp = mean(Pamp(validI,1:2), 2, 'omitnan');
shortAmp = mean(Pamp(validI,3:4), 2, 'omitnan');
fprintf(fid, '\nLate harmonic long/short anisotropy:\n');
fprintf(fid, '  mean(|P_f0|_long/|P_f0|_short)=%.8g\n', mean(longAmp./shortAmp,'omitnan'));
fprintf(fid, '  mean(|P_f0|_long-|P_f0|_short)=%.8g\n', mean(longAmp-shortAmp,'omitnan'));
fprintf(fid, '\nRepresentative stations:\n');
for ii = iStations
    fprintf(fid, '  i=%3d X=%.8g  Long+ rms=%.8g amp=%.8g phase=%.8g | Short+ rms=%.8g amp=%.8g phase=%.8g\n', ...
        ii, Xline(ii,1), Prms(ii,1), Pamp(ii,1), Pphase(ii,1), Prms(ii,3), Pamp(ii,3), Pphase(ii,3));
end
fprintf(fid, '\nLate-window three-subwindow RMS convergence, rows=windows, columns Long+/Short+ at stations:\n');
for w = 1:size(conv.rms,1)
    fprintf(fid, '  window %d:', w);
    for s = 1:numel(iStations)
        fprintf(fid, ' i=%d L=%.8g S=%.8g;', iStations(s), conv.rms(w,s,1), conv.rms(w,s,2));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\nWall maps: max(P_rms)=%.8g, max(|P_f0|)=%.8g\n', max(PwallRms(:)), max(PwallAmp(:)));
fclose(fid);
end

function write_late_csv(outDir, Xline, labels, Pmean, Prms, Pamp, Pphase)
fid = fopen(fullfile(outDir, 'late_axis_harmonic_profiles.csv'), 'w');
fprintf(fid, 'i,x,label,P_mean_late,P_rms_late,P_harmonic_amp,P_harmonic_phase_rad\n');
for kk = 1:numel(labels)
    for i = 1:size(Prms,1)
        fprintf(fid, '%d,%.10g,%s,%.10g,%.10g,%.10g,%.10g\n', ...
            i, Xline(i,kk), labels{kk}, Pmean(i,kk), Prms(i,kk), Pamp(i,kk), Pphase(i,kk));
    end
end
fclose(fid);
end

function plot_axis_quantity(outDir, Xline, labels, A, ttl, ylab, fname)
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

function plot_axis_phase(outDir, Xline, labels, phase, f0)
fig = figure('Visible','off','Position',[100,100,1100,720]);
for kk = [1 3]
    plot(Xline(:,kk), unwrap(phase(:,kk)), 'LineWidth', 2.2); hold on;
end
xlabel('X'); ylabel('phase (rad)');
title(sprintf('Late-window harmonic phase of P'', f0 = %.4g', f0));
legend({labels{1}, labels{3}}, 'Location','best');
grid on; box on; style_pub(gca);
print(fig, fullfile(outDir, 'late_axis_pressure_harmonic_phase.png'), '-dpng', '-r220');
close(fig);
end

function plot_wall_map(outDir, Xwall, A, kList, labels, ttl, cbLabel, fname)
xLine = mean(Xwall,2,'omitnan');
plotA = A;
plotCbLabel = cbLabel;
if contains(fname, 'harmonic_amplitude')
    plotA = A * 1e3;
    plotCbLabel = '|P''_{f0}| (10^{-3})';
end
fig = figure('Visible','off','Position',[100,100,1900,1180]);
imagesc(xLine, 1:size(plotA,2), plotA');
axis xy; hold on;
styles = {'-k','-k','--k','--k'};
for kk = 1:numel(kList)
    plot(xLine, kList(kk)*ones(size(xLine)), styles{kk}, 'LineWidth', 3.4);
    labelY = min(max(kList(kk)+8, 14), size(A,2)-10);
    text(max(xLine)*0.982, labelY, labels{kk}, 'HorizontalAlignment','right', ...
        'FontName','Times New Roman','FontSize',38,'FontWeight','bold', ...
        'BackgroundColor','w','EdgeColor','k','Margin',8,'LineWidth',1.0);
end
xlabel('X'); ylabel('circumferential index k'); title(ttl);
cb = colorbar; ylabel(cb, plotCbLabel);
grid on; box on; style_pub(gca);
set(gca, 'FontSize',38, 'LineWidth',1.8);
set(get(gca,'XLabel'), 'FontSize',44);
set(get(gca,'YLabel'), 'FontSize',44);
set(get(gca,'Title'), 'FontSize',38);
set(cb, 'FontName','Times New Roman','FontSize',38,'FontWeight','bold','LineWidth',1.8);
set(get(cb,'Label'), 'FontName','Times New Roman','FontSize',44,'FontWeight','bold');
print(fig, fullfile(outDir, fname), '-dpng', '-r220');
close(fig);
end

function plot_wall_phase(outDir, Xwall, phase, kList, labels, f0)
xLine = mean(Xwall,2,'omitnan');
fig = figure('Visible','off','Position',[100,100,1900,1180]);
imagesc(xLine, 1:size(phase,2), phase');
axis xy; colormap(hsv); caxis([-pi pi]); hold on;
styles = {'-k','-k','--k','--k'};
for kk = 1:numel(kList)
    plot(xLine, kList(kk)*ones(size(xLine)), styles{kk}, 'LineWidth', 3.4);
    labelY = min(max(kList(kk)+8, 14), size(phase,2)-10);
    text(max(xLine)*0.982, labelY, labels{kk}, 'HorizontalAlignment','right', ...
        'FontName','Times New Roman','FontSize',38,'FontWeight','bold', ...
        'BackgroundColor','w','EdgeColor','k','Margin',8,'LineWidth',1.0);
end
xlabel('X'); ylabel('circumferential index k');
title(sprintf('Late-window harmonic phase of wall pressure perturbation, f0=%.4g', f0));
cb = colorbar; ylabel(cb, 'phase (rad)');
grid on; box on; style_pub(gca);
set(gca, 'FontSize',38, 'LineWidth',1.8);
set(get(gca,'XLabel'), 'FontSize',44);
set(get(gca,'YLabel'), 'FontSize',44);
set(get(gca,'Title'), 'FontSize',38);
set(cb, 'FontName','Times New Roman','FontSize',38,'FontWeight','bold','LineWidth',1.8);
set(get(cb,'Label'), 'FontName','Times New Roman','FontSize',44,'FontWeight','bold');
print(fig, fullfile(outDir, 'late_wall_pressure_harmonic_phase_map.png'), '-dpng', '-r220');
close(fig);
end

function plot_late_time_fit(outDir, t, Pline, Xline, labels, iStations, lateIds, f0)
fig = figure('Visible','off','Position',[100,100,1300,850]);
for s = 1:numel(iStations)
    ii = iStations(s);
    subplot(numel(iStations),1,s);
    tt = t(lateIds);
    yLong = squeeze(Pline(lateIds,ii,1));
    yShort = squeeze(Pline(lateIds,ii,3));
    [ts, yLongS] = smooth_time_series(tt, yLong);
    [~, yShortS] = smooth_time_series(tt, yShort);
    plot(ts, yLongS, '-', 'LineWidth', 2.2); hold on;
    plot(ts, yShortS, '--', 'LineWidth', 2.2);
    plot(tt, yLong, 'o', 'MarkerSize', 4, 'LineWidth', 1.0, 'HandleVisibility','off');
    plot(tt, yShort, 's', 'MarkerSize', 4, 'LineWidth', 1.0, 'HandleVisibility','off');
    ylabel('P''');
    title(sprintf('Late-window P'' time history, i=%d, X=%.4g', ii, Xline(ii,1)));
    legend({labels{1}, labels{3}}, 'Location','best');
    grid on; box on; style_pub(gca);
end
xlabel('time');
print(fig, fullfile(outDir, 'late_pressure_perturbation_time_histories.png'), '-dpng', '-r220');
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
ts = linspace(min(t0), max(t0), max(300, numel(t0)*14))';
ys = interp1(t0, y0, ts, 'pchip');
end

function plot_convergence(outDir, conv, iStations)
fig = figure('Visible','off','Position',[100,100,1100,700]);
for s = 1:numel(iStations)
    subplot(1,numel(iStations),s);
    plot(1:size(conv.rms,1), squeeze(conv.rms(:,s,1)), '-o', 'LineWidth',2.0); hold on;
    plot(1:size(conv.rms,1), squeeze(conv.rms(:,s,2)), '--s', 'LineWidth',2.0);
    xlabel('late sub-window'); ylabel('P'' RMS');
    title(sprintf('i=%d', iStations(s)));
    legend(conv.labels, 'Location','best');
    grid on; box on; style_pub(gca);
end
print(fig, fullfile(outDir, 'late_subwindow_rms_convergence.png'), '-dpng', '-r220');
close(fig);
end

function style_pub(ax)
set(ax, 'FontName','Times New Roman', 'FontSize',22, 'FontWeight','bold', 'LineWidth',1.3);
set(get(ax,'XLabel'), 'FontName','Times New Roman', 'FontSize',26, 'FontWeight','bold');
set(get(ax,'YLabel'), 'FontName','Times New Roman', 'FontSize',26, 'FontWeight','bold');
set(get(ax,'Title'), 'FontName','Times New Roman', 'FontSize',26, 'FontWeight','bold');
lgd = findobj(ancestor(ax,'figure'), 'Type','Legend');
for i = 1:numel(lgd)
    set(lgd(i), 'FontName','Times New Roman', 'FontSize',14, 'FontWeight','bold');
end
end
