function results = interpolateFlowFieldFolder(dataFolder, filePattern, gridStruct, colStruct, options)
% interpolateFlowFieldFolder
%
% Interpolates scattered CSV flow-field data onto a structured grid.
%
% INPUTS:
%   dataFolder   - folder containing CSV files
%   filePattern  - filename pattern (e.g. '*.csv')
%   gridStruct   - struct with:
%         .xmin, .xmax, .ymin, .ymax
%         .dx, .dy
%   colStruct    - struct defining column indices:
%         .x, .y (required)
%         .u, .v (optional)
%         .alpha (optional)
%   options      - struct with optional fields:
%         .maxFiles
%         .interpMethod   (default 'natural')
%         .extrapMethod   (default 'nearest')
%
% OUTPUT:
%   results struct containing:
%         .fields  (interpolated data)
%         .xgrid, .ygrid
%         .files
%

%% ----------------------------
% Defaults
%% ----------------------------
if nargin < 5
    options = struct;
end

if ~isfield(options, 'maxFiles')
    options.maxFiles = inf;
end

if ~isfield(options, 'interpMethod')
    options.interpMethod = 'natural';
end

if ~isfield(options, 'extrapMethod')
    options.extrapMethod = 'nearest';
end

%% ----------------------------
% Create structured grid
%% ----------------------------
xgrid = gridStruct.xmin : gridStruct.dx : gridStruct.xmax;
ygrid = gridStruct.ymin : gridStruct.dy : gridStruct.ymax;

[X,Y] = meshgrid(xgrid, ygrid);

Nx = numel(xgrid);
Ny = numel(ygrid);

%% ----------------------------
% Get file list
%% ----------------------------
files = dir(fullfile(dataFolder, filePattern));
nFiles = min(numel(files), options.maxFiles);

if nFiles == 0
    error('No matching CSV files found.');
end

fprintf('Grid size: %d x %d\n', Nx, Ny);
fprintf('Files to process: %d\n\n', nFiles);

%% ----------------------------
% Determine which variables exist
%% ----------------------------
hasU     = isfield(colStruct,'u');
hasV     = isfield(colStruct,'v');
hasAlpha = isfield(colStruct,'alpha');

%% ----------------------------
% Preallocate dynamically
%% ----------------------------
fields = struct;

if hasU
    fields.u = nan(nFiles, Ny, Nx);
end

if hasV
    fields.v = nan(nFiles, Ny, Nx);
end

if hasAlpha
    fields.alpha = nan(nFiles, Ny, Nx);
end

%% ============================
% MAIN LOOP
%% ============================
for it = 1:nFiles

    fprintf('Processing %d / %d\n', it, nFiles);

    filename = fullfile(dataFolder, files(it).name);

    data = readmatrix(filename);

    % Extract coordinates
    x = data(:, colStruct.x);
    y = data(:, colStruct.y);

    % Restrict domain before interpolation (faster)
    mask = x >= gridStruct.xmin & x <= gridStruct.xmax & ...
           y >= gridStruct.ymin & y <= gridStruct.ymax;

    x = x(mask);
    y = y(mask);

    % ---- U velocity
    if hasU
        u = data(mask, colStruct.u);
        F = scatteredInterpolant(x,y,u, ...
                                 options.interpMethod, ...
                                 options.extrapMethod);
        fields.u(it,:,:) = F(X,Y);
    end

    % ---- V velocity
    if hasV
        v = data(mask, colStruct.v);
        F = scatteredInterpolant(x,y,v, ...
                                 options.interpMethod, ...
                                 options.extrapMethod);
        fields.v(it,:,:) = F(X,Y);
    end

    % ---- Alpha
    if hasAlpha
        alpha = data(mask, colStruct.alpha);
        F = scatteredInterpolant(x,y,alpha, ...
                                 options.interpMethod, ...
                                 options.extrapMethod);
        fields.alpha(it,:,:) = F(X,Y);
    end

end

%% ----------------------------
% Output struct
%% ----------------------------
results.fields = fields;
results.xgrid  = xgrid;
results.ygrid  = ygrid;
results.files  = {files(1:nFiles).name};

fprintf('\nInterpolation complete.\n');

end