function [slip_hr, res_hr, valid_idx] = upscaleSSLIP(grain_ebsdID, ebsd_grain)
% upscaleSSLIP mathematically inverts the spatial coarse-graining of SSLIP results.
% 
% This function extracts the spatial downsampling metadata (sps, In, Im) inherently 
% embedded in grain_ebsdID by preprocessSSLIP, mathematically upscales the raw matrices 
% using repelem, embeds them into the uncropped bounds, and flattens them for easy 
% injection into the global arrays.
%
% Outputs:
%   slip_hr   - The upscaled, high-res flattened slipIDcor matrix [N_pixels x N_sys]
%   res_hr    - The upscaled, high-res flattened residualEeff array [N_pixels x 1]
%   valid_idx - The exact 1D linear memory indices of the valid pixels in the global EBSD map

    % 1. Extract the perfect 1-to-1 mapping arrays precomputed by preprocessSSLIP.m
    id_hr = grain_ebsdID.prop.upscale_id_hr;
    idx_lr = grain_ebsdID.prop.upscale_idx_lr;
    
    % 2. Extract the solver results safely as pure MATLAB arrays (bypassing MTEX truncation bugs)
    slip_raw = grain_ebsdID.prop.slipIDcor;
    if ndims(slip_raw) == 3
        N_sys = size(slip_raw, 3);
        slip_lr = reshape(slip_raw, [], N_sys);
    else
        slip_lr = slip_raw;
    end
    res_lr = grain_ebsdID.prop.residualEeff(:);
    
    % 3. Map directly! The precomputed tables inherently upsample and crop-align the results perfectly!
    slip_hr = slip_lr(idx_lr, :);
    res_hr = res_lr(idx_lr);
    valid_idx = id_hr;
end
