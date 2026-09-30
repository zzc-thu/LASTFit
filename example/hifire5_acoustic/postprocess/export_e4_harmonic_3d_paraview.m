function export_e4_harmonic_3d_paraview()
% Export late-periodic f0 harmonic pressure amplitude and phase on the 3-D grid.

rawRoot = getenv('LASTFIT_HIFIRE_RAW_DIR');
if isempty(rawRoot)
    error(['Set LASTFIT_HIFIRE_RAW_DIR to the directory containing ', ...
        'E-2, E-3, E-4, and RESU (for example, Un_Hifire_fine).']);
end
caseDir = fullfile(rawRoot, 'E-4');
pertDir = fullfile(caseDir, 'Pert');
outDir = fullfile(caseDir, 'paraview_harmonic_f0_3d');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

cfg = read_config(fullfile(caseDir, 'Config.cfg'));
dtStep = get_cfg_number(cfg, 'FixedDeltaTime', 1e-7);
stepsWrite = get_cfg_number(cfg, 'StepsWriteData', 2000);
dtSave = dtStep * stepsWrite;

masters = dir(fullfile(pertDir, 'Pert_Result_*.pvts'));
if isempty(masters)
    masters = dir(fullfile(pertDir, 'Result_*.pvts'));
end
if isempty(masters)
    error('No Pert master PVTS files found in %s.', pertDir);
end
[frameIds, order] = sort(extract_frame_ids({masters.name}));
masters = masters(order);
nFrames = numel(masters);
timeAll = (0:nFrames-1) * dtSave;

lateStart = ceil(nFrames/2);
lateIdx = lateStart:nFrames;
timeLate = timeAll(lateIdx);

% The dominant peak found from the full-window temporal spectrum.
f0 = 666.6666666667;
omega = 2*pi*f0;
phaseVector = exp(-1i * omega * timeLate(:));

firstMaster = fullfile(pertDir, masters(1).name);
masterInfo = parse_pvts(firstMaster);
nPieces = numel(masterInfo.pieces);

fprintf('Pert frames: %d, late-periodic frames: %d..%d (%d frames)\n', ...
    nFrames, lateStart, nFrames, numel(lateIdx));
fprintf('dt_save = %.12g, f0 = %.8g\n', dtSave, f0);
fprintf('PVTS pieces: %d\n', nPieces);

outPieceNames = cell(nPieces, 1);
outExtents = cell(nPieces, 1);
wallPieceNames = {};
wallExtents = {};
wallContribs = {};
globalAmpMax = 0;
wallAmpMax = 0;

for ip = 1:nPieces
    fprintf('Processing piece %d / %d...\n', ip, nPieces);
    y = [];
    points = [];
    pieceExtent = '';

    for it = 1:numel(lateIdx)
        mf = fullfile(pertDir, masters(lateIdx(it)).name);
        info = parse_pvts(mf);
        pieceRel = info.pieces(ip).source;
        piecePath = fullfile(pertDir, pieceRel);
        piece = read_vts_piece(piecePath, 'P_pert');

        if it == 1
            points = piece.points;
            pieceExtent = piece.extent;
            y = zeros(numel(piece.values), numel(lateIdx));
        end
        y(:, it) = piece.values(:);
    end

    y = y - mean(y, 2);
    coeff = (2 / numel(lateIdx)) * (y * phaseVector);
    amp = abs(coeff);
    phaseRad = angle(coeff);
    phaseDeg = phaseRad * 180 / pi;
    realPart = real(coeff);
    imagPart = imag(coeff);
    globalAmpMax = max(globalAmpMax, max(amp));

    outPieceNames{ip} = sprintf('Harmonic_f0_666p67_Part_%05d.vts', ip-1);
    outExtents{ip} = pieceExtent;
    write_vts_piece(fullfile(outDir, outPieceNames{ip}), pieceExtent, points, ...
        amp, phaseRad, phaseDeg, realPart, imagPart);

    [hasWall, wallExtent, wallIdx] = wall_slice_indices(pieceExtent, 1);
    if hasWall
        wallPieceNames{end+1, 1} = sprintf('Harmonic_f0_666p67_Wall_Part_%05d.vts', ip-1); %#ok<AGROW>
        wallExtents{end+1, 1} = wallExtent; %#ok<AGROW>
        write_vts_piece(fullfile(outDir, wallPieceNames{end}), wallExtent, points(wallIdx, :), ...
            amp(wallIdx), phaseRad(wallIdx), phaseDeg(wallIdx), realPart(wallIdx), imagPart(wallIdx));
        wallContribs{end+1, 1} = make_wall_contrib(wallExtent, points(wallIdx, :), ...
            realPart(wallIdx), imagPart(wallIdx)); %#ok<AGROW>
        wallAmpMax = max(wallAmpMax, max(amp(wallIdx)));
    end
end

write_pvts_master(fullfile(outDir, 'Harmonic_f0_666p67_3D.pvts'), ...
    masterInfo.wholeExtent, outExtents, outPieceNames);
write_pvts_master(fullfile(outDir, 'Harmonic_f0_666p67_Wall_3D.pvts'), ...
    wall_whole_extent(masterInfo.wholeExtent, 1), wallExtents, wallPieceNames);
write_smooth_wall_outputs(outDir, wallContribs, f0);
write_paraview_style_script(outDir);
write_correspondence_note(outDir, f0);

write_summary(fullfile(outDir, 'harmonic_f0_3d_summary.txt'), ...
    nFrames, lateStart, nFrames, dtSave, f0, numel(lateIdx), globalAmpMax, wallAmpMax, frameIds);

fprintf('Done. Open in ParaView:\n%s\n', fullfile(outDir, 'Harmonic_f0_666p67_3D.pvts'));
fprintf('Wall surface file:\n%s\n', fullfile(outDir, 'Harmonic_f0_666p67_Wall_3D.pvts'));
fprintf('Smoothed single-wall file:\n%s\n', fullfile(outDir, 'Harmonic_f0_666p67_Wall_Smooth.vts'));
end

function cfg = read_config(pathName)
cfg = containers.Map();
if ~exist(pathName, 'file')
    return;
end
txt = fileread(pathName);
lines = regexp(txt, '\r?\n', 'split');
for i = 1:numel(lines)
    line = strtrim(regexprep(lines{i}, '[#!%].*$', ''));
    if isempty(line)
        continue;
    end
    tok = regexp(line, '^\s*([^=\s]+)\s*=\s*([^,\s]+)', 'tokens', 'once');
    if ~isempty(tok)
        cfg(tok{1}) = tok{2};
    end
end
end

function value = get_cfg_number(cfg, key, defaultValue)
if isKey(cfg, key)
    value = str2double(cfg(key));
    if ~isnan(value)
        return;
    end
end
value = defaultValue;
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

function piece = read_vts_piece(pathName, varName)
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
    extentTok = regexp(txt, 'WholeExtent="([^"]+)"', 'tokens', 'once');
end
if isempty(extentTok)
    error('Cannot find piece extent in %s.', pathName);
end
pieceExtent = extentTok{1};
nPoints = count_points_from_extent(pieceExtent);

pointsBlock = regexp(txt, '<Points>\s*(.*?)\s*</Points>', 'tokens', 'once');
if isempty(pointsBlock)
    error('Cannot find Points block in %s.', pathName);
end
pointsArray = regexp(pointsBlock{1}, '<DataArray([^>]*)/?>', 'tokens', 'once');
if isempty(pointsArray)
    error('Cannot find Points DataArray in %s.', pathName);
end
points = read_data_array(raw, txt, markerPos, marker, pointsArray{1}, pathName, 3*nPoints);
if mod(numel(points), 3) ~= 0
    error('Point coordinate count is not divisible by 3 in %s.', pathName);
end
points = reshape(points, 3, []).';

valueExpr = ['<DataArray([^>]*)Name="' regexptranslate('escape', varName) '"([^>]*)/?>'];
valueParts = regexp(txt, valueExpr, 'tokens', 'once');
if isempty(valueParts)
    error('Cannot find variable %s in %s.', varName, pathName);
end
values = read_data_array(raw, txt, markerPos, marker, [valueParts{1} ' ' valueParts{2}], pathName, nPoints);

if size(points, 1) ~= numel(values)
    error('Point/value size mismatch in %s: %d points, %d values.', ...
        pathName, size(points, 1), numel(values));
end

piece.extent = pieceExtent;
piece.points = points;
piece.values = values;
end

function nPoints = count_points_from_extent(extentText)
nums = sscanf(extentText, '%f');
if numel(nums) ~= 6
    error('Invalid VTK extent: %s.', extentText);
end
nPoints = prod(nums(2:2:6) - nums(1:2:5) + 1);
end

function [hasWall, wallExtent, idx] = wall_slice_indices(extentText, wallJ)
nums = sscanf(extentText, '%f');
if numel(nums) ~= 6
    error('Invalid VTK extent: %s.', extentText);
end
iMin = nums(1); iMax = nums(2);
jMin = nums(3); jMax = nums(4);
kMin = nums(5); kMax = nums(6);
hasWall = wallJ >= jMin && wallJ <= jMax;
wallExtent = '';
idx = [];
if ~hasWall
    return;
end
nI = iMax - iMin + 1;
nJ = jMax - jMin + 1;
nK = kMax - kMin + 1;
jLocal = wallJ - jMin + 1;
idx = zeros(nI*nK, 1);
c = 0;
for kk = 1:nK
    for ii = 1:nI
        c = c + 1;
        idx(c) = ii + (jLocal-1)*nI + (kk-1)*nI*nJ;
    end
end
wallExtent = sprintf('%d %d %d %d %d %d', iMin, iMax, wallJ, wallJ, kMin, kMax);
end

function extentText = wall_whole_extent(wholeExtent, wallJ)
nums = sscanf(wholeExtent, '%f');
if numel(nums) ~= 6
    error('Invalid whole extent: %s.', wholeExtent);
end
nums(3) = wallJ;
nums(4) = wallJ;
extentText = sprintf('%d %d %d %d %d %d', nums(1), nums(2), nums(3), nums(4), nums(5), nums(6));
end

function contrib = make_wall_contrib(wallExtent, points, realPart, imagPart)
nums = sscanf(wallExtent, '%f');
nI = nums(2) - nums(1) + 1;
nK = nums(6) - nums(5) + 1;
contrib.extent = nums(:).';
contrib.points = points;
contrib.realPart = realPart(:);
contrib.imagPart = imagPart(:);
contrib.nI = nI;
contrib.nK = nK;
end

function write_smooth_wall_outputs(outDir, wallContribs, f0)
if isempty(wallContribs)
    warning('No wall contributions found; skipping smooth wall output.');
    return;
end

allExt = zeros(numel(wallContribs), 6);
for c = 1:numel(wallContribs)
    allExt(c, :) = wallContribs{c}.extent;
end
iMin = min(allExt(:, 1)); iMax = max(allExt(:, 2));
kMin = min(allExt(:, 5)); kMax = max(allExt(:, 6));
nI = iMax - iMin + 1;
nK = kMax - kMin + 1;

X = zeros(nI, nK);
Y = zeros(nI, nK);
Z = zeros(nI, nK);
R = zeros(nI, nK);
Q = zeros(nI, nK);
C = zeros(nI, nK);

for c = 1:numel(wallContribs)
    item = wallContribs{c};
    e = item.extent;
    ii = (e(1):e(2)) - iMin + 1;
    kk = (e(5):e(6)) - kMin + 1;
    ptsX = reshape(item.points(:, 1), item.nI, item.nK);
    ptsY = reshape(item.points(:, 2), item.nI, item.nK);
    ptsZ = reshape(item.points(:, 3), item.nI, item.nK);
    valsR = reshape(item.realPart, item.nI, item.nK);
    valsQ = reshape(item.imagPart, item.nI, item.nK);
    X(ii, kk) = X(ii, kk) + ptsX;
    Y(ii, kk) = Y(ii, kk) + ptsY;
    Z(ii, kk) = Z(ii, kk) + ptsZ;
    R(ii, kk) = R(ii, kk) + valsR;
    Q(ii, kk) = Q(ii, kk) + valsQ;
    C(ii, kk) = C(ii, kk) + 1;
end

mask = C > 0;
X(mask) = X(mask) ./ C(mask);
Y(mask) = Y(mask) ./ C(mask);
Z(mask) = Z(mask) ./ C(mask);
R(mask) = R(mask) ./ C(mask);
Q(mask) = Q(mask) ./ C(mask);
tmp = nearest_fill(X, mask); X(~mask) = tmp(~mask);
tmp = nearest_fill(Y, mask); Y(~mask) = tmp(~mask);
tmp = nearest_fill(Z, mask); Z(~mask) = tmp(~mask);
tmp = nearest_fill(R, mask); R(~mask) = tmp(~mask);
tmp = nearest_fill(Q, mask); Q(~mask) = tmp(~mask);

complexRaw = R + 1i * Q;
ampRaw = abs(complexRaw);
phaseRaw = angle(complexRaw);

smoothPasses = 4;
complexSmooth = smooth_complex_wall(complexRaw, smoothPasses);
ampSmooth = abs(complexSmooth);
phaseSmooth = angle(complexSmooth);
phaseUnwrappedX = unwrap(phaseSmooth, [], 1);
phaseContinuous = continuous_phase_display(phaseUnwrappedX);

points = zeros(nI*nK, 3);
arrays = {
    'P_f0_amplitude', ampRaw(:);
    'P_f0_phase_rad', phaseRaw(:);
    'P_f0_phase_deg', phaseRaw(:) * 180 / pi;
    'P_f0_real', real(complexRaw(:));
    'P_f0_imag', imag(complexRaw(:));
    'P_f0_amplitude_smooth', ampSmooth(:);
    'P_f0_amplitude_x1e3_smooth', ampSmooth(:) * 1e3;
    'P_f0_phase_rad_smooth', phaseSmooth(:);
    'P_f0_phase_deg_smooth', phaseSmooth(:) * 180 / pi;
    'P_f0_phase_unwrapped_x_rad_smooth', phaseUnwrappedX(:);
    'P_f0_phase_continuous_display', phaseContinuous(:);
    'P_f0_phase_continuous_display_deg', phaseContinuous(:) * 180 / pi;
    'P_f0_real_smooth', real(complexSmooth(:));
    'P_f0_imag_smooth', imag(complexSmooth(:));
    'wall_i_index', repmat((iMin:iMax).', nK, 1);
    'wall_k_index', reshape(repmat(kMin:kMax, nI, 1), [], 1);
    };
idx = 0;
for kk = 1:nK
    for ii = 1:nI
        idx = idx + 1;
        points(idx, :) = [X(ii, kk), Y(ii, kk), Z(ii, kk)];
    end
end

extent = sprintf('%d %d 1 1 %d %d', iMin, iMax, kMin, kMax);
write_vts_piece_arrays(fullfile(outDir, 'Harmonic_f0_666p67_Wall_Smooth.vts'), ...
    extent, points, arrays, 'P_f0_amplitude_smooth');

summaryPath = fullfile(outDir, 'harmonic_f0_wall_smooth_summary.txt');
fid = fopen(summaryPath, 'w');
if fid >= 0
    fprintf(fid, 'Smoothed single-wall f0 harmonic output\n');
    fprintf(fid, 'f0: %.12g\n', f0);
    fprintf(fid, 'wall extent: %s\n', extent);
    fprintf(fid, 'grid size: nI=%d, nK=%d\n', nI, nK);
    fprintf(fid, 'smoothing: complex coefficient 3x3 weighted passes = %d\n', smoothPasses);
    fprintf(fid, 'raw wall max amplitude: %.12g\n', max(ampRaw(:)));
    fprintf(fid, 'smooth wall max amplitude: %.12g\n', max(ampSmooth(:)));
    fprintf(fid, '\nRecommended ParaView arrays:\n');
    fprintf(fid, 'Amplitude: P_f0_amplitude_smooth, use range [0, %.4g]\n', max(ampSmooth(:)));
    fprintf(fid, 'Wrapped phase: P_f0_phase_rad_smooth, use range [-pi, pi]\n');
    fprintf(fid, 'Unwrapped streamwise phase: P_f0_phase_unwrapped_x_rad_smooth\n');
    fprintf(fid, 'Best visual phase: P_f0_phase_continuous_display, use symmetric data range\n');
    fclose(fid);
end
end

function out = nearest_fill(A, mask)
out = zeros(size(A));
if any(mask(:))
    out(:) = mean(A(mask), 'omitnan');
end
end

function Zs = smooth_complex_wall(Z, passes)
Zs = Z;
for p = 1:passes
    Zs = smooth_once_periodic_k(Zs);
end
end

function phaseDisplay = continuous_phase_display(phaseUnwrappedX)
phaseDisplay = phaseUnwrappedX;
phaseDisplay = phaseDisplay - mean(phaseDisplay(:), 'omitnan');
for p = 1:3
    phaseDisplay = smooth_once_periodic_k(phaseDisplay);
end
lim = prctile(abs(phaseDisplay(:)), 98);
if isfinite(lim) && lim > 0
    phaseDisplay = max(min(phaseDisplay, lim), -lim);
end
end

function B = smooth_once_periodic_k(A)
kernelI = [1 4 6 4 1];
kernelK = [1 4 6 4 1];
B = zeros(size(A));
weightSum = sum(kernelI) * sum(kernelK);
for di = -2:2
    Ai = shift_i_clamped(A, di);
    for dk = -2:2
        B = B + kernelI(di+3) * kernelK(dk+3) * circshift(Ai, [0, dk]);
    end
end
B = B / weightSum;
end

function B = shift_i_clamped(A, di)
B = zeros(size(A));
if di < 0
    B(1:end+di, :) = A(1-di:end, :);
    B(end+di+1:end, :) = repmat(A(end, :), -di, 1);
elseif di > 0
    B(1+di:end, :) = A(1:end-di, :);
    B(1:di, :) = repmat(A(1, :), di, 1);
else
    B = A;
end
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

function write_vts_piece(pathName, extent, points, amp, phaseRad, phaseDeg, realPart, imagPart)
fid = fopen(pathName, 'w');
if fid < 0
    error('Cannot open %s for writing.', pathName);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

arrays = {
    'P_f0_amplitude', single(amp(:));
    'P_f0_phase_rad', single(phaseRad(:));
    'P_f0_phase_deg', single(phaseDeg(:));
    'P_f0_real', single(realPart(:));
    'P_f0_imag', single(imagPart(:));
    };
pointData = single(reshape(points.', [], 1));
offsets = zeros(size(arrays, 1) + 1, 1, 'uint64');
runningOffset = uint64(0);
for i = 1:numel(offsets)
    offsets(i) = runningOffset;
    if i == 1
        nBytes = numel(pointData) * 4;
    else
        nBytes = numel(arrays{i-1, 2}) * 4;
    end
    runningOffset = runningOffset + uint64(8 + nBytes);
end

fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n');
fprintf(fid, '  <StructuredGrid WholeExtent="%s">\n', extent);
fprintf(fid, '    <Piece Extent="%s">\n', extent);
fprintf(fid, '      <Points>\n');
fprintf(fid, '        <DataArray type="Float32" NumberOfComponents="3" format="appended" offset="%d"/>\n', offsets(1));
fprintf(fid, '      </Points>\n');
fprintf(fid, '      <PointData Scalars="P_f0_amplitude">\n');
for i = 1:size(arrays, 1)
    fprintf(fid, '        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="appended" offset="%d"/>\n', ...
        arrays{i, 1}, offsets(i+1));
end
fprintf(fid, '      </PointData>\n');
fprintf(fid, '      <CellData>\n');
fprintf(fid, '      </CellData>\n');
fprintf(fid, '    </Piece>\n');
fprintf(fid, '  </StructuredGrid>\n');
fprintf(fid, '  <AppendedData encoding="raw">_');
write_appended_block(fid, pointData);
for i = 1:size(arrays, 1)
    write_appended_block(fid, arrays{i, 2});
end
fprintf(fid, '\n  </AppendedData>\n');
fprintf(fid, '</VTKFile>\n');
end

function write_vts_piece_arrays(pathName, extent, points, arrays, scalarName)
fid = fopen(pathName, 'w');
if fid < 0
    error('Cannot open %s for writing.', pathName);
end
cleanup = onCleanup(@() fclose(fid));

pointData = single(reshape(points.', [], 1));
offsets = zeros(size(arrays, 1) + 1, 1, 'uint64');
runningOffset = uint64(0);
offsets(1) = runningOffset;
runningOffset = runningOffset + uint64(8 + numel(pointData) * 4);
for i = 1:size(arrays, 1)
    offsets(i+1) = runningOffset;
    runningOffset = runningOffset + uint64(8 + numel(arrays{i, 2}) * 4);
end

fprintf(fid, '<VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n');
fprintf(fid, '  <StructuredGrid WholeExtent="%s">\n', extent);
fprintf(fid, '    <Piece Extent="%s">\n', extent);
fprintf(fid, '      <Points>\n');
fprintf(fid, '        <DataArray type="Float32" NumberOfComponents="3" format="appended" offset="%d"/>\n', offsets(1));
fprintf(fid, '      </Points>\n');
fprintf(fid, '      <PointData Scalars="%s">\n', scalarName);
for i = 1:size(arrays, 1)
    fprintf(fid, '        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="appended" offset="%d"/>\n', ...
        arrays{i, 1}, offsets(i+1));
end
fprintf(fid, '      </PointData>\n');
fprintf(fid, '      <CellData>\n');
fprintf(fid, '      </CellData>\n');
fprintf(fid, '    </Piece>\n');
fprintf(fid, '  </StructuredGrid>\n');
fprintf(fid, '  <AppendedData encoding="raw">_');
write_appended_block(fid, pointData);
for i = 1:size(arrays, 1)
    write_appended_block(fid, single(arrays{i, 2}(:)));
end
fprintf(fid, '\n  </AppendedData>\n');
fprintf(fid, '</VTKFile>\n');
end

function write_appended_block(fid, data)
nBytes = uint64(numel(data) * 4);
fwrite(fid, nBytes, 'uint64');
fwrite(fid, data, 'single');
end

function write_pvts_master(pathName, wholeExtent, extents, sources)
fid = fopen(pathName, 'w');
if fid < 0
    error('Cannot open %s for writing.', pathName);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="PStructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n');
fprintf(fid, '  <PStructuredGrid WholeExtent="%s" GhostLevel="0">\n', wholeExtent);
fprintf(fid, '    <PPointData Scalars="P_f0_amplitude">\n');
fprintf(fid, '      <PDataArray type="Float32" Name="P_f0_amplitude" NumberOfComponents="1"/>\n');
fprintf(fid, '      <PDataArray type="Float32" Name="P_f0_phase_rad" NumberOfComponents="1"/>\n');
fprintf(fid, '      <PDataArray type="Float32" Name="P_f0_phase_deg" NumberOfComponents="1"/>\n');
fprintf(fid, '      <PDataArray type="Float32" Name="P_f0_real" NumberOfComponents="1"/>\n');
fprintf(fid, '      <PDataArray type="Float32" Name="P_f0_imag" NumberOfComponents="1"/>\n');
fprintf(fid, '    </PPointData>\n');
fprintf(fid, '    <PCellData>\n');
fprintf(fid, '    </PCellData>\n');
fprintf(fid, '    <PPoints>\n');
fprintf(fid, '      <PDataArray type="Float32" NumberOfComponents="3"/>\n');
fprintf(fid, '    </PPoints>\n');
for i = 1:numel(sources)
    fprintf(fid, '    <Piece Extent="%s" Source="%s"/>\n', extents{i}, sources{i});
end
fprintf(fid, '  </PStructuredGrid>\n');
fprintf(fid, '</VTKFile>\n');
end

function write_paraview_style_script(outDir)
pathName = fullfile(outDir, 'load_harmonic_wall_smooth_style.py');
fid = fopen(pathName, 'w');
if fid < 0
    warning('Cannot write ParaView style script.');
    return;
end
cleanup = onCleanup(@() fclose(fid));
wallFile = strrep(fullfile(outDir, 'Harmonic_f0_666p67_Wall_Smooth.vts'), '\', '/');
fprintf(fid, 'from paraview.simple import *\n');
fprintf(fid, 'paraview.simple._DisableFirstRenderCameraReset()\n\n');
fprintf(fid, 'src = OpenDataFile(r"%s")\n', wallFile);
fprintf(fid, 'view = GetActiveViewOrCreate("RenderView")\n');
fprintf(fid, 'disp = Show(src, view)\n');
fprintf(fid, 'disp.Representation = "Surface"\n');
fprintf(fid, 'ColorBy(disp, ("POINTS", "P_f0_amplitude_smooth"))\n');
fprintf(fid, 'ampLUT = GetColorTransferFunction("P_f0_amplitude_smooth")\n');
fprintf(fid, 'ampLUT.RGBPoints = [0.0, 0.10, 0.20, 0.70, 0.0015, 0.55, 0.75, 1.00, 0.0030, 1.00, 0.93, 0.60, 0.0045, 0.70, 0.05, 0.05]\n');
fprintf(fid, 'ampLUT.ColorSpace = "RGB"\n');
fprintf(fid, 'ampLUT.NanColor = [0.6, 0.6, 0.6]\n');
fprintf(fid, 'disp.RescaleTransferFunctionToDataRange(False, True)\n');
fprintf(fid, 'ampLUT.RescaleTransferFunction(0.0, 0.0045)\n');
fprintf(fid, 'ampBar = GetScalarBar(ampLUT, view)\n');
fprintf(fid, 'ampBar.Title = "P_f0 amplitude (smooth)"\n');
fprintf(fid, 'ampBar.ComponentTitle = ""\n');
fprintf(fid, 'ampBar.TitleFontSize = 16\n');
fprintf(fid, 'ampBar.LabelFontSize = 14\n');
fprintf(fid, 'view.Background = [1.0, 1.0, 1.0]\n');
fprintf(fid, 'view.OrientationAxesVisibility = 1\n');
fprintf(fid, 'view.Update()\n');
fprintf(fid, 'ResetCamera(view)\n\n');
fprintf(fid, '# To display phase with the same smoothed wall file, run these lines in ParaView Python Shell:\n');
fprintf(fid, '# ColorBy(disp, ("POINTS", "P_f0_phase_continuous_display"))\n');
fprintf(fid, '# phaseLUT = GetColorTransferFunction("P_f0_phase_continuous_display")\n');
fprintf(fid, '# phaseLUT.RGBPoints = [-6.0, 0.12, 0.22, 0.85, 0.0, 0.95, 0.95, 0.95, 6.0, 0.75, 0.05, 0.05]\n');
fprintf(fid, '# phaseLUT.ColorSpace = "Diverging"\n');
fprintf(fid, '# phaseLUT.RescaleTransferFunction(-6.0, 6.0)\n');
fprintf(fid, '# GetScalarBar(phaseLUT, view).Title = "P_f0 phase trend (continuous)"\n');
end

function write_correspondence_note(outDir, f0)
pathName = fullfile(outDir, 'paraview_wall_correspondence.md');
fid = fopen(pathName, 'w', 'n', 'UTF-8');
if fid < 0
    warning('Cannot write ParaView correspondence note.');
    return;
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# ParaView 壁面谐波场与二维傅里叶图的对应关系\n\n');
fprintf(fid, '本目录中的 `Harmonic_f0_666p67_Wall_Smooth.vts` 是将后半段压力扰动 `P_pert` 在 `f0 = %.6g Hz` 上的 Fourier 复系数重新贴回三维壁面 `j=1` 后得到的单一壁面结构网格。\n\n', f0);
fprintf(fid, '二维图 `late_wall_pressure_harmonic_amplitude_map.png` 的横坐标 `X` 与该 VTK 壁面上的空间坐标 `X` 对应，纵坐标 circumferential index `k` 与 VTK 数组 `wall_k_index` 对应。因此二维图中的任意点 `(X,k)`，在 ParaView 中就是壁面上相同 `X` 位置和相同 `wall_k_index` 的点。\n\n');
fprintf(fid, '二维幅值图中的 `|p''_{f0}|` 对应 ParaView 数组 `P_f0_amplitude_smooth`。若需要和二维图色标的 `10^{-3}` 标注一致，也可以显示 `P_f0_amplitude_x1e3_smooth`。\n\n');
fprintf(fid, '二维相位图中的 phase(rad) 对应 ParaView 数组 `P_f0_phase_rad_smooth`。该数组仍然限制在 `[-pi, pi]`，适合与二维相位图直接对比，但在三维表面上会保留相位包裹跳变。若主要目的是展示连续传播趋势，建议显示 `P_f0_phase_continuous_display`；若需要完整未包裹相位，可显示 `P_f0_phase_unwrapped_x_rad_smooth`。\n\n');
fprintf(fid, '四条参考轴线的周向索引为：Long axis +: `k=41`，Long axis -: `k=121`，Short axis +: `k=1`，Short axis -: `k=81`。在 ParaView 中可以用 `Threshold` 或 `Plot Over Line` 结合 `wall_k_index` 提取这些轴线，与二维图中的黑色实线/虚线标注对应。\n\n');
fprintf(fid, '本次新增的平滑不是直接平滑相位角，而是先对 Fourier 复系数 `P_f0_real + i P_f0_imag` 在壁面 `(i,k)` 上做 5x5 加权平滑，再重新计算幅值和相位，因此可以避免 `-pi/pi` 相位跳变处的错误平均。\n');
end

function write_summary(pathName, nFrames, lateStart, lateEnd, dtSave, f0, nLate, globalAmpMax, wallAmpMax, frameIds)
fid = fopen(pathName, 'w');
if fid < 0
    error('Cannot open %s for writing.', pathName);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '3-D f0 harmonic export for E-4 Pert field\n');
fprintf(fid, 'Total Pert frames: %d\n', nFrames);
fprintf(fid, 'Frame ids: %d .. %d\n', frameIds(1), frameIds(end));
fprintf(fid, 'Late periodic frame indices used in time order: %d .. %d\n', lateStart, lateEnd);
fprintf(fid, 'Late periodic frame count: %d\n', nLate);
fprintf(fid, 'Save interval dt_save: %.12g\n', dtSave);
fprintf(fid, 'Sampling frequency: %.12g\n', 1/dtSave);
fprintf(fid, 'Exported dominant frequency f0: %.12g\n', f0);
fprintf(fid, 'Maximum P_f0_amplitude in exported 3-D field: %.12g\n', globalAmpMax);
fprintf(fid, 'Maximum P_f0_amplitude on exported wall surface: %.12g\n', wallAmpMax);
fprintf(fid, '\nArrays written to ParaView:\n');
fprintf(fid, 'P_f0_amplitude: pressure perturbation harmonic amplitude at f0\n');
fprintf(fid, 'P_f0_phase_rad: pressure perturbation harmonic phase in radians\n');
fprintf(fid, 'P_f0_phase_deg: pressure perturbation harmonic phase in degrees\n');
fprintf(fid, 'P_f0_real/P_f0_imag: complex harmonic coefficient components\n');
end
