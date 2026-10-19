function [prepData, ebsdID] = preprocessSSLIP(ebsd, DeformationData, cfg_prep)
% preprocessSSLIP Handles data filtering, coarse-graining, and gradient extraction 
% before running the SSLIP optimizer.

if nargin < 3, cfg_prep = struct; end

X = ebsd.x;
Y = ebsd.y;

hasGradients = isfield(DeformationData, 'Hxx') && isfield(DeformationData, 'Hxy') && isfield(DeformationData, 'Hyx') && isfield(DeformationData, 'Hyy');
hasDisps = isfield(DeformationData, 'U') && isfield(DeformationData, 'V');

if ~hasGradients && ~hasDisps
    error('DeformationData must contain either U and V, or Hxx, Hxy, Hyx, Hyy');
end

% Initialize dummy structure for spatial coordinates
crs = struct();

% --- PROCESS DISPLACEMENTS ---
if hasDisps
    fprintf('Calculating displacement gradients from provided displacements...\n');
    U = DeformationData.U;
    V = DeformationData.V;
    
    if cfg_prep.filterSize ~= 0
        filt_opts.filt_std = cfg_prep.filterSize;
        U(U == 0) = NaN;
        V(V == 0) = NaN;
        data = filterDisplacements(U, V, filt_opts);
    else
        data.U = U;
        data.V = V;
    end
    
    crs = coarsegrainDisp(data.U, data.V, X(1,:), Y(:,1)', cfg_prep.coarsegrain);
    
    if ~hasGradients
        % Calculate gradients numerically
        [Hxx, Hxy] = gradient(crs.f, crs.pixelsize(1), crs.pixelsize(2));
        [Hyx, Hyy] = gradient(crs.g, crs.pixelsize(1), crs.pixelsize(2));
    end
    
    prepData.U = crs.f;
    prepData.V = crs.g;
end

% --- PROCESS GRADIENTS ---
if hasGradients
    fprintf('Using provided displacement gradients...\n');
    Hxx = DeformationData.Hxx; Hxy = DeformationData.Hxy;
    Hyx = DeformationData.Hyx; Hyy = DeformationData.Hyy;
    
    % If raw grid gradients were passed (same size as original EBSD)
    if isequal(size(Hxx), size(ebsd.x))
        if cfg_prep.filterSize ~= 0
            filt_opts.filt_std = cfg_prep.filterSize;
            fH1 = filterDisplacements(Hxx, Hxy, filt_opts);
            fH2 = filterDisplacements(Hyx, Hyy, filt_opts);
            Hxx = fH1.U; Hxy = fH1.V;
            Hyx = fH2.U; Hyy = fH2.V;
        end
        if cfg_prep.coarsegrain > 1
            crsH1 = coarsegrainDisp(Hxx, Hxy, X(1,:), Y(:,1)', cfg_prep.coarsegrain);
            crsH2 = coarsegrainDisp(Hyx, Hyy, X(1,:), Y(:,1)', cfg_prep.coarsegrain);
            Hxx = crsH1.f; Hxy = crsH1.g;
            Hyx = crsH2.f; Hyy = crsH2.g;
        end
    end
    
    if ~hasDisps
        % Need spatial grid info even if Displacements aren't provided
        if cfg_prep.coarsegrain > 1
            crs = coarsegrainDisp(zeros(size(ebsd.x)), zeros(size(ebsd.y)), X(1,:), Y(:,1)', cfg_prep.coarsegrain);
        else
            crs.X = X; crs.Y = Y;
        end
        prepData.U = zeros(size(Hxx));
        prepData.V = zeros(size(Hxx));
    end
end

% Create EBSD mapping object
ebsdID = dummyEBSDSimple(ebsd.orientations(1), crs.X, crs.Y);

prepData.Hxx = Hxx; prepData.Hxy = Hxy;
prepData.Hyx = Hyx; prepData.Hyy = Hyy;
prepData.Eeff = calcEffectiveE(Hxx, Hxy, Hyx, Hyy);

% Store preprocessed data in EBSD object
ebsdID.prop.U = prepData.U;
ebsdID.prop.V = prepData.V;
ebsdID.prop.Hxx = Hxx; ebsdID.prop.Hxy = Hxy;
ebsdID.prop.Hyx = Hyx; ebsdID.prop.Hyy = Hyy;
ebsdID.prop.Eeff = prepData.Eeff;

% -------------------------------------------------------------------------
% Precompute the EXACT 1-to-1 mathematical upscaling mapping!
% By computing this here, we completely isolate all downstream scripts 
% from ANY spatial math, cropping logic, or grid reshaping!
% -------------------------------------------------------------------------
if isfield(crs, 'sps') && crs.sps > 1
    In = crs.In_bounds;
    Im = crs.Im_bounds;
    sps = crs.sps;
    
    % The original exact MTEX .id values of the un-cropped region (using global_id if stashed)
    if isfield(ebsd.prop, 'global_id')
        id_kept = ebsd.prop.global_id(In, Im);
    else
        id_kept = ebsd.id(In, Im);
    end
    
    % The grid size of the low-res solver
    n_lr = length(In) / sps;
    m_lr = length(Im) / sps;
    
    % A 2D grid representing the native 1D linear memory indices of the low-res solver
    idx_lr_grid = reshape(1:(n_lr * m_lr), n_lr, m_lr);
    
    % Upscale the index grid perfectly using the superpixel size
    idx_upscaled = repelem(idx_lr_grid, sps, sps);
    
    % Filter out invalid/NaN IDs just to be safe
    valid = ~isnan(id_kept(:)) & (id_kept(:) > 0);
    
    % Attach the perfectly aligned 1D lookup tables to the object!
    ebsdID.prop.upscale_id_hr = id_kept(valid);
    ebsdID.prop.upscale_idx_lr = idx_upscaled(valid);
else
    % If no coarsegraining was applied, the mapping is just 1-to-1 identity
    if isfield(ebsd.prop, 'global_id')
        id_full = ebsd.prop.global_id;
    else
        id_full = ebsd.id;
    end
    valid = ~isnan(id_full(:)) & (id_full(:) > 0);
    ebsdID.prop.upscale_id_hr = id_full(valid);
    ebsdID.prop.upscale_idx_lr = find(valid);
end

end
