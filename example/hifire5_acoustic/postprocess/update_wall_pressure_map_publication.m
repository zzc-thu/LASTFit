function update_wall_pressure_map_publication()
% Export a publication-style wall pressure map to normal_profiles_separate.

scriptDir = fileparts(mfilename('fullpath'));
caseReleaseDir = fileparts(scriptDir);
rootDir = getenv('LASTFIT_HIFIRE_RAW_DIR');
if isempty(rootDir)
    error(['Set LASTFIT_HIFIRE_RAW_DIR to the directory containing ', ...
        'E-2, E-3, E-4, and RESU (for example, Un_Hifire_fine).']);
end
resuDir = fullfile(rootDir, 'RESU');
outDir = fullfile(caseReleaseDir, 'output', 'manuscript_fine', ...
    'base_flow_wall_pressure');
pvtsFile = fullfile(resuDir, 'Result_00000179.pvts');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

info = parse_pvts_local(pvtsFile);
dims = info.whole(2:2:6) - info.whole(1:2:5) + 1;
ni = dims(1);
nk = dims(3);
wallJGlobal = 1;

Xwall = nan(ni, nk);
Pwall = nan(ni, nk);

for p = 1:numel(info.pieces)
    piece = info.pieces(p);
    e = piece.extent;
    if wallJGlobal < e(3) || wallJGlobal > e(4)
        continue;
    end

    data = read_vts_piece_wall(fullfile(resuDir, piece.source), {'P'});
    ii = (e(1):e(2)) - info.whole(1) + 1;
    kk = (e(5):e(6)) - info.whole(5) + 1;
    jj = wallJGlobal - e(3) + 1;
    Xwall(ii, kk) = double(squeeze(data.X(:, jj, :)));
    Pwall(ii, kk) = double(squeeze(data.P(:, jj, :)));
end

xLine = mean(Xwall, 2, 'omitnan');
kLine = 1:nk;

longPos = 41;
longNeg = 121;
shortPos = 1;
shortNeg = 81;

fig = figure('Visible','off','Position',[100,100,1320,900]);
imagesc(xLine, kLine, Pwall');
axis xy;
hold on;
plot_axis_line(longPos, xLine, 'Long axis +', '-k');
plot_axis_line(longNeg, xLine, 'Long axis -', '-k');
plot_axis_line(shortPos, xLine, 'Short axis +', '--k');
plot_axis_line(shortNeg, xLine, 'Short axis -', '--k');

xlabel('X');
ylabel('Circumferential index k');
title('Wall pressure');
cb = colorbar;
ylabel(cb, 'P');
grid on;
box on;

tickSize = 16;
labelSize = 32;
set(gca, 'FontName', 'Times New Roman', 'FontSize', tickSize, 'FontWeight', 'bold', ...
    'LineWidth', 1.4, 'Box', 'on');
set(get(gca, 'XLabel'), 'FontName', 'Times New Roman', 'FontSize', labelSize, 'FontWeight', 'bold');
set(get(gca, 'YLabel'), 'FontName', 'Times New Roman', 'FontSize', labelSize, 'FontWeight', 'bold');
set(get(gca, 'Title'), 'FontName', 'Times New Roman', 'FontSize', labelSize, 'FontWeight', 'bold');
set(cb, 'FontName', 'Times New Roman', 'FontSize', tickSize, 'FontWeight', 'bold', 'LineWidth', 1.2);
set(get(cb, 'YLabel'), 'FontName', 'Times New Roman', 'FontSize', labelSize, 'FontWeight', 'bold');

print(fig, fullfile(outDir, 'wall_pressure_map.png'), '-dpng', '-r300');
close(fig);
fprintf('Updated %s\n', fullfile(outDir, 'wall_pressure_map.png'));
end

function plot_axis_line(k, xLine, labelText, styleSpec)
plot(xLine, k * ones(size(xLine)), styleSpec, 'LineWidth', 1.8);
text(max(xLine) * 0.985, k + 2.5, labelText, ...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
    'FontName','Times New Roman', 'FontSize',24, 'FontWeight','bold', ...
    'Color','k', 'BackgroundColor','w', 'Margin',1);
end

function info = parse_pvts_local(pvtsFile)
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

function data = read_vts_piece_wall(fn, names)
fid = fopen(fn, 'r', 'ieee-le');
if fid < 0
    error('Cannot open %s', fn);
end
bytes = fread(fid, inf, 'uint8=>uint8');
fclose(fid);
pat = uint8('<AppendedData encoding="raw">_');
idx = strfind(bytes', pat);
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

for i = 1:numel(names)
    nm = names{i};
    expr = ['<DataArray[^>]*Name="' nm '"[^>]*offset="\s*(\d+)'];
    tok = regexp(header, expr, 'tokens', 'once');
    off = str2double(tok{1});
    start = base + off + 8;
    vals = typecast(bytes(start:start + npts*4 - 1), 'single');
    data.(nm) = reshape(vals, ni, nj, nk);
end
end
