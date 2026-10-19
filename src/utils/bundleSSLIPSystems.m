function [bundles, metric_vals, leftover_idx] = bundleSSLIPSystems(ebsdID, sSLocal, opt)
% bundleSSLIPSystems Groups slip systems based on user-defined similarity criteria.
%
% Inputs:
%   ebsdID   - EBSD object containing .prop.slipIDcor
%   sSLocal  - List of local or reference slip systems
%   opt      - Options struct governing the bundling logic:
%              .activityPercentile (default 85)
%              .activityThreshold  (default 0.01) - Systems below this threshold are thrown into leftover_idx and not bundled. If set to 0, all systems are bundled.
%              .bundleKinematic    (default true) - Match 2D deformation gradient. Evaluated in specimen coordinates (2D projection).
%              .bundleCommonSlipSystems (default true) - Match exact n and b Miller indices. Evaluated in crystal coordinates (requires opt.ori to map specimen vectors back to crystal frame).
%              .bundleBurgers      (default false) - Match Burgers vector angle. Evaluated directly on sSLocal.b (specimen coordinates).
%              .bundleNormalTrace  (default false) - Match slip trace angle on the observation plane. Evaluated directly on sSLocal.trace (specimen coordinates).
%              .bundleNormal       (default false) - Match slip plane normal angle. Evaluated directly on sSLocal.n (specimen coordinates).
%              .bundleStrainThresh (default 0.2)
%              .bundleAngleThresh  (default 1e-4 radians)
%              .bundleTraceAngleThresh (default 3 degrees) - Tolerance for trace matching (used for option bundleNormalTrace)
%              .ignoreZeroActivity (default true) - Strip out near-zero noise pixels (gamma <= 1e-5) before calculating activity percentiles. If false, zeroes are included and percentiles may drop to 0.
% Outputs:
%   bundles      - Cell array of grouped system indices
%   metric_vals  - Computed activity metrics for all systems
%   leftover_idx - Indices of systems that did not meet activity threshold

if nargin < 3, opt = struct; end

% --- DEFAULTS ---
if ~isfield(opt, 'activityPercentile'), opt.activityPercentile = 85; end
if ~isfield(opt, 'activityThreshold'),  opt.activityThreshold = 0; end
if ~isfield(opt, 'bundleStrainThresh'), opt.bundleStrainThresh = 0.2; end
if ~isfield(opt, 'bundleAngleThresh'),  opt.bundleAngleThresh = 1e-4; end
if ~isfield(opt, 'bundleTraceAngleThresh'),  opt.bundleTraceAngleThresh = 5*degree; end

% Toggles
if ~isfield(opt, 'bundleKinematic'),         opt.bundleKinematic = true; end
if ~isfield(opt, 'bundleCommonSlipSystems'), opt.bundleCommonSlipSystems = true; end
if ~isfield(opt, 'bundleBurgers'),           opt.bundleBurgers = false; end
if ~isfield(opt, 'bundleNormalTrace'),       opt.bundleNormalTrace = false; end
if ~isfield(opt, 'bundleNormal'),            opt.bundleNormal = false; end
if ~isfield(opt, 'ignoreZeroActivity'),      opt.ignoreZeroActivity = false; end

% --- PRECOMPUTE COMPARISON MATRICES ---
num_total_sys = length(sSLocal);
is_crystallographic = (1:num_total_sys <= num_total_sys); % Update this if pseudo-systems are appended

if opt.bundleKinematic
    Hslip = sSLocal.deformationTensor.matrix;
    if isfield(opt, 'enableRotation') && opt.enableRotation
        Hslip(:,:,1,end+1) = [0 -1 0; 1 0 0; 0 0 0];
        if isfield(opt, 'posConstr') && opt.posConstr
            Hslip(:,:,1,end+1) = [0 1 0; -1 0 0; 0 0 0];
        end
        is_crystallographic = [true(1, num_total_sys), false(1, size(Hslip,4)-num_total_sys)];
        num_total_sys = size(Hslip, 4);
    end
    Hslip11 = reshape(Hslip(1,1,:), 1, []);
    Hslip12 = reshape(Hslip(1,2,:), 1, []);
    Hslip21 = reshape(Hslip(2,1,:), 1, []);
    Hslip22 = reshape(Hslip(2,2,:), 1, []);
    A_mat = [Hslip11; Hslip12; Hslip21; Hslip22]; 
    if isfield(opt, 'normalizeInplane') && opt.normalizeInplane
        A_mat = A_mat ./ sqrt(sum(A_mat.^2, 1));
    end
else
    if isfield(opt, 'enableRotation') && opt.enableRotation
        num_sys_actual = num_total_sys + 1;
        if isfield(opt, 'posConstr') && opt.posConstr
            num_sys_actual = num_sys_actual + 1;
        end
        is_crystallographic = [true(1, num_total_sys), false(1, num_sys_actual-num_total_sys)];
        num_total_sys = num_sys_actual;
    end
end

% Let MTEX isolate the indexed EBSD pixels
ebsd_indexed = ebsdID('indexed');

% Map the local .id of the indexed pixels back to the global EBSD object's .id array.
% This returns the exact 1D row indices corresponding to the active pixels, 
% completely bypassing any row-major / column-major logical mask flattening bugs.
[~, valid_idx] = ismember(ebsd_indexed.id, ebsdID.id);

% Filter out any invalid mappings (just in case)
valid_idx = valid_idx(valid_idx > 0);

% Extract the subset using robust 1D linear row indexing
gamma_abs = abs(ebsdID.prop.slipIDcor(valid_idx, :));

num_sys_data = size(gamma_abs, 2);
metric_vals = zeros(1, num_sys_data);
max_vals = zeros(1, num_sys_data);
maxGamma = zeros(1, num_sys_data);

for i = 1:num_sys_data
    gamma_i = gamma_abs(:, i);
    
    if opt.ignoreZeroActivity
        valid_g = gamma_i(~isnan(gamma_i) & ~isinf(gamma_i) & gamma_i > 1e-5);
    else
        valid_g = gamma_i(~isnan(gamma_i) & ~isinf(gamma_i));
    end
    
    if isempty(valid_g)
        metric_vals(i) = 0;
        max_vals(i) = 0;
    else
        metric_vals(i) = prctile(valid_g, opt.activityPercentile);
        max_vals(i) = max(valid_g);
    end
    maxGamma(i) = max(gamma_i);
end

% If metric_vals are uniformly crushed to zero, sort by max_vals instead to preserve order!
if max(metric_vals) == 0 && max(max_vals) > 0
    [~, sort_sys_idx] = sort(max_vals, 'descend');
else
    [~, sort_sys_idx] = sort(metric_vals, 'descend');
end

% --- BUNDLING LOGIC ---
bundles = {};
rep_sys_list = [];
leftover_idx = [];

for i = 1:length(sort_sys_idx)
    sys_i = sort_sys_idx(i);
    
    if metric_vals(sys_i) < opt.activityThreshold
        leftover_idx = [leftover_idx, sys_i];
        continue;
    end
    
    matched_bundle = 0;
    
    if ~isempty(rep_sys_list)
        % Can only do crystallographic checks on true slip systems
        is_cryst_i = is_crystallographic(sys_i);
        
        for b_idx = 1:length(rep_sys_list)
            rep = rep_sys_list(b_idx);
            is_cryst_rep = is_crystallographic(rep);
            
            match_kinematic = true;
            match_cryst     = true;
            match_burgers   = true;
            match_trace     = true;
            match_normal    = true;
            
            % 1. Kinematic Check
            if opt.bundleKinematic
                if norm(A_mat(:, rep) - A_mat(:, sys_i)) > opt.bundleStrainThresh
                    match_kinematic = false;
                end
            end
            
            % If it's a valid crystallographic system, do crystallographic checks
            if is_cryst_i && is_cryst_rep
                % 2. Strict Crystallographic Check (Needs opt.ori)
                if opt.bundleCommonSlipSystems && isfield(opt, 'ori') && ~isempty(opt.ori)
                    miller_rep_n = round(inv(opt.ori) * sSLocal(rep).n);
                    miller_sys_n = round(inv(opt.ori) * sSLocal(sys_i).n);
                    miller_rep_b = round(inv(opt.ori) * sSLocal(rep).b, 'uvw');
                    miller_sys_b = round(inv(opt.ori) * sSLocal(sys_i).b, 'uvw');
                    
                    if angle(miller_rep_n, miller_sys_n, 'noSymmetry', 'antipodal') > 1e-5 || ...
                       angle(miller_rep_b, miller_sys_b, 'noSymmetry', 'antipodal') > 1e-5
                        match_cryst = false;
                    end
                elseif opt.bundleCommonSlipSystems && (~isfield(opt, 'ori') || isempty(opt.ori))
                    warning('bundleSSLIPSystems:MissingOri', 'bundleCommonSlipSystems is true, but opt.ori is missing. Cannot check system equivalence. Omitted.');
                end
                
                % 3. Burgers Vector Check
                if opt.bundleBurgers
                    b_ang = angle(sSLocal(rep).b, sSLocal(sys_i).b, 'noSymmetry', 'antipodal');
                    if b_ang > opt.bundleAngleThresh
                        match_burgers = false;
                    end
                end
                
                % 4. Trace Angle Check
                if opt.bundleNormalTrace
                    t_ang = angle(sSLocal(rep).trace, sSLocal(sys_i).trace);
                    if t_ang > opt.bundleTraceAngleThresh
                        match_trace = false;
                    end
                end
                
                % 5. Normal Angle Check
                if opt.bundleNormal
                    n_ang = angle(sSLocal(rep).n, sSLocal(sys_i).n, 'noSymmetry', 'antipodal');
                    if n_ang > opt.bundleAngleThresh
                        match_normal = false;
                    end
                end
            else
                % If one is pseudo-slip and the other isn't, they never match via crystallography.
                if is_cryst_i ~= is_cryst_rep
                    match_cryst = false; match_burgers = false; match_trace = false; match_normal = false;
                else
                    % If BOTH are pseudo-slip (e.g. rotation), they only match if they are literally the same system
                    % which doesn't happen across unique indices.
                    if sys_i ~= rep
                        match_cryst = false;
                    end
                end
            end
            
            % % DEBUGGER FOR SYSTEM n AND m
            % n = 11;
            % m = 47
            % if (sys_i == n && rep == m) || (sys_i == m && rep == n)
            %     fprintf('\n--- DEBUG BUNDLE: Sys %d vs Sys %d ---\n', rep, sys_i);
            %     if opt.bundleKinematic
            %         fprintf('Kinematic Match: %d (Norm diff: %f vs Thresh: %f)\n', match_kinematic, norm(A_mat(:, rep) - A_mat(:, sys_i)), opt.bundleStrainThresh);
            %     end
            %     if is_cryst_i && is_cryst_rep
            %         if opt.bundleCommonSlipSystems && isfield(opt, 'ori') && ~isempty(opt.ori)
            %             fprintf('Common Systems : %d (N ang: %f deg, B ang: %f deg)\n', match_cryst, ...
            %                 angle(miller_rep_n, miller_sys_n, 'noSymmetry', 'antipodal')/degree, ...
            %                 angle(miller_rep_b, miller_sys_b, 'noSymmetry', 'antipodal')/degree);
            %         end
            %         if opt.bundleBurgers
            %             fprintf('Burgers Match  : %d (Ang: %f deg vs Thresh: %f deg)\n', match_burgers, angle(sSLocal(rep).b, sSLocal(sys_i).b, 'noSymmetry', 'antipodal')/degree, opt.bundleAngleThresh/degree);
            %         end
            %         if opt.bundleNormalTrace
            %             fprintf('Trace Match    : %d (Ang: %f deg vs Thresh: %f deg)\n', match_trace, angle(sSLocal(rep).trace, sSLocal(sys_i).trace)/degree, opt.bundleTraceAngleThresh/degree);
            %         end
            %         if opt.bundleNormal
            %             fprintf('Normal Match   : %d (Ang: %f deg vs Thresh: %f deg)\n', match_normal, angle(sSLocal(rep).n, sSLocal(sys_i).n, 'noSymmetry', 'antipodal')/degree, opt.bundleAngleThresh/degree);
            %         end
            %     end
            %     fprintf('---------------------------------------\n');
            % end
            
            % Combine flags (Logical AND logic as specified)
            if match_kinematic && match_cryst && match_burgers && match_trace && match_normal
                matched_bundle = b_idx;
                break;
            end
        end
    end
    
    if matched_bundle > 0
        bundles{matched_bundle} = [bundles{matched_bundle}, sys_i];
    else
        bundles{end+1} = sys_i;
        rep_sys_list = [rep_sys_list, sys_i];
    end
end
end
