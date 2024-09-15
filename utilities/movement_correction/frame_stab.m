function m = frame_stab(m)
% Frame stabilization using separable Gaussian filtering
%   Jinghao Lu, 10/28/2018
%
% This function applies spatial and temporal Gaussian filtering to stabilize
% the frames in the 'reg' variable of the matfile 'm'. It processes the data
% in batches to manage memory usage and avoid overloading the system.

    %%% Initialize variables %%%
    [pixh, pixw, nf] = size(m, 'reg'); % Get dimensions of the 'reg' variable
    sigmas = [1, 1, 1]; % Standard deviations for the Gaussian kernels in each dimension
    hsizes = 2 * ceil(2 * sigmas) + 1; % Kernel sizes based on the sigmas
    [hcol, hrow, hslc] = createSeparableGaussianKernel(sigmas, hsizes); % Create Gaussian kernels
    
    % Determine the data type and compute an estimate of data size
    ttype = class(m.reg(1, 1, 1));
    stype = parse_type(ttype); % Get size in bytes of the data type
    nsize = pixh * pixw * nf * stype * 10; % Heuristic size for batch computation
    
    %%% First pass: Spatial filtering along columns and rows %%%
    nbatch = batch_compute(nsize); % Compute initial number of batches based on data size
    nbatch = min(nbatch, nf); % Ensure nbatch does not exceed the number of frames (nf)
    idbatch = unique(round(linspace(0, nf, nbatch + 1))); % Indices to split frames into batches
    nbatch = length(idbatch) - 1; % Update the actual number of batches
    
    for i = 1: nbatch
        % Define start and end indices for the batch
        idx_start = idbatch(i) + 1;
        idx_end = idbatch(i + 1);
        if idx_start > idx_end
            continue; % Skip invalid ranges where start index exceeds end index
        end
        % Extract the batch of frames
        tmp = m.reg(1: pixh, 1: pixw, idx_start: idx_end);
        % Convolve with the column and row kernels (spatial filtering)
        tmp = convn(tmp, hcol, 'same');
        tmp = convn(tmp, hrow, 'same');
        % Write the filtered frames back into 'm.reg'
        m.reg(1: pixh, 1: pixw, idx_start: idx_end) = tmp;
    end
    
    %%% Second pass: Temporal filtering along frames %%%
    nbatch = batch_compute(nsize); % Recompute number of batches for temporal filtering
    nbatch = min(nbatch, pixh); % Ensure nbatch does not exceed the number of rows (pixh)
    idbatch = unique(round(linspace(0, pixh, nbatch + 1))); % Indices to split rows into batches
    nbatch = length(idbatch) - 1; % Update the actual number of batches
    
    for i = 1: nbatch
        % Define start and end indices for the batch
        idx_start = idbatch(i) + 1;
        idx_end = idbatch(i + 1);
        if idx_start > idx_end
            continue; % Skip invalid ranges where start index exceeds end index
        end
        % Extract the batch of rows across all frames (spatial dimension)
        tmp = m.reg(idx_start: idx_end, 1: pixw, 1: nf);
        % Convolve with the slice (temporal) kernel (temporal filtering)
        tmp = convn(tmp, hslc, 'same');
        % Write the filtered data back into 'm.reg'
        m.reg(idx_start: idx_end, 1: pixw, 1: nf) = tmp;
    end
end

function [hcol, hrow, hslc] = createSeparableGaussianKernel(sigma, hsize)
% Create separable Gaussian kernels for convolution in each dimension.
% Inputs:
%   sigma - Standard deviations for the Gaussian kernels [sigma_x, sigma_y, sigma_z]
%   hsize - Sizes of the kernels in each dimension [size_x, size_y, size_z]
% Outputs:
%   hcol - Column kernel for the x-dimension (columns)
%   hrow - Row kernel for the y-dimension (rows)
%   hslc - Slice kernel for the z-dimension (slices or frames)

    % Check if the kernel is isotropic (same sigma and hsize in all dimensions)
    isIsotropic = all(sigma == sigma(1)) && all(hsize == hsize(1));

    % Create the column kernel (1D Gaussian for x-dimension)
    hcol = images.internal.createGaussianKernel(sigma(1), hsize(1));

    if isIsotropic
        % If isotropic, reuse the same kernel for all dimensions
        hrow = hcol;
        hslc = hcol;
    else
        % Otherwise, create separate kernels for each dimension
        hrow = images.internal.createGaussianKernel(sigma(2), hsize(2));
        hslc = images.internal.createGaussianKernel(sigma(3), hsize(3));
    end

    % Reshape the kernels to be suitable for convolution
    hrow = reshape(hrow, 1, hsize(2)); % Reshape for y-dimension convolution
    hslc = reshape(hslc, 1, 1, hsize(3)); % Reshape for z-dimension convolution
end



% ORIGINAL BELOW
% function m = frame_stab(m)
% % frame stabilization 
% %   Jinghao Lu, 10/28/2018
% 
%     %%% initialize %%%
%     [pixh, pixw, nf] = size(m, 'reg');
%     sigmas = [1, 1, 1];
%     hsizes = 2 * ceil(2 * sigmas) + 1;
%     [hcol,hrow,hslc] = createSeparableGaussianKernel(sigmas, hsizes);
%     ttype = class(m.reg(1, 1, 1));
%     stype = parse_type(ttype);
%     nsize = pixh * pixw * nf * stype * 10; %%% heuristic size of algorithm %%%
% 
%     %%% first spatial %%%
%     nbatch = batch_compute(nsize);
%     idbatch = round(linspace(0, nf, nbatch + 1));
% 
%     for i = 1: nbatch
%         tmp = m.reg(1: pixh, 1: pixw, idbatch(i) + 1: idbatch(i + 1));
%         tmp = convn(tmp, hcol, 'same');
%         tmp = convn(tmp, hrow, 'same');
%         m.reg(1: pixh, 1: pixw, idbatch(i) + 1: idbatch(i + 1)) = tmp;
%     end
% 
%     %%% second temporal %%%
%     nbatch = batch_compute(nsize);
%     idbatch = round(linspace(0, pixh, nbatch + 1));
% 
%     for i = 1: nbatch
%         tmp = m.reg(idbatch(i) + 1: idbatch(i + 1), 1: pixw, 1: nf);
%         tmp = convn(tmp, hslc, 'same');
%         m.reg(idbatch(i) + 1: idbatch(i + 1), 1: pixw, 1: nf) = tmp;
%     end
% end
% 
% function [hcol,hrow,hslc] = createSeparableGaussianKernel(sigma, hsize)
% 
% isIsotropic = all(sigma==sigma(1)) && all(hsize==hsize(1));
% 
% hcol = images.internal.createGaussianKernel(sigma(1), hsize(1));
% 
% if isIsotropic
%     hrow = hcol;
%     hslc = hcol;
% else
%     hrow = images.internal.createGaussianKernel(sigma(2), hsize(2));
%     hslc = images.internal.createGaussianKernel(sigma(3), hsize(3));
% end
% 
% hrow = reshape(hrow, 1, hsize(2));
% hslc = reshape(hslc, 1, 1, hsize(3));
% 
% end
