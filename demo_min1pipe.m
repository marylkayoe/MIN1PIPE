%%% demo of the full MIN1PIPE %%%

% Clear variables from the workspace to ensure a clean environment
clear all;

%% Session-Specific Parameter Initialization %%
Fsi = 60;             % Original sampling rate (Hz)
Fsi_new = 60;         % New sampling rate after temporal downsampling (Hz)
spatialr = 0.5;         % Spatial downsampling factor (1 means no downsampling)
se = 9;               % Structure element size for background removal
ismc = 0;             % Motion correction flag (1 for applying motion correction, 0 otherwise)
flag = 2;             % Seed selection method (1 for automatic, 2 for manual)
% isvis = true;       % Visualization flag (unused here)
ifpost = false;       % Post-processing flag (set to true if post-processing is desired)

%% Main Program %%
% Run the main pipeline function 'min1pipe' with the specified parameters
% This function processes the imaging data and saves the results
[fname, frawname, fregname] = min1pipe(Fsi, Fsi_new, spatialr, se, ismc, flag);

%% Plot Some Images %%
% Load the results from 'min1pipe'
load(fname)  % This loads variables like 'imaxn', 'imaxy', 'imax', 'sigfn', 'roifn', 'seedsfn', etc.

figure(1)
clf

%%% Plot Raw Maximum Projection %%%
subplot(2, 3, 1, 'align')
imagesc(imaxn)  % 'imaxn' is the maximum projection of the raw data
axis square
title('Raw')

%%% Plot Neural Enhanced Image Before Motion Correction %%%
subplot(2, 3, 2, 'align')
imagesc(imaxy)  % 'imaxy' is the maximum projection after initial enhancement but before motion correction
axis square
title('Before MC')

%%% Plot Neural Enhanced Image After Motion Correction %%%
subplot(2, 3, 3, 'align')
imagesc(imax)   % 'imax' is the maximum projection after motion correction
axis square
title('After MC')

%%% Plot Contours of Detected Cells %%%
subplot(2, 3, 4, 'align')
if ~isempty(sigfn) && ~isempty(roifn) && ~isempty(seedsfn)
    % If cells were detected, plot their contours
    plot_contour(roifn, sigfn, seedsfn, imax, pixh, pixw)
else
    % If no cells were detected, display a message
    axis off
    title('No Cells Found')
end
axis square

%%% Plot All Identified Traces %%%
subplot(2, 3, [5 6])
if ~isempty(sigfn)
    % Normalize and plot the fluorescence traces of detected cells
    sigt = sigfn;
    for i = 1: size(sigt, 1)
        sigt(i, :) = normalize(sigt(i, :));
    end
    plot((sigt + (1: size(sigt, 1))')')
    axis tight
    title('Traces')
else
    % If no traces are available, display a message
    axis off
    title('No Traces to Display')
end

%% Making a TIFF Stack %%
% Load the raw and registered data files
load(fname)          % Ensure variables like 'pixh' and 'pixw' are loaded
mraw = matfile(frawname);
mreg = matfile(fregname);

% Determine the data type of the frames based on motion correction flag
if ismc
    % If motion correction was applied, use 'frame_all' (motion-corrected data)
    ttype = class(mraw.frame_all(1, 1, 1));
else
    % If motion correction was not applied, use 'frame_allt' (preprocessed data without motion correction)
    ttype = class(mraw.frame_allt(1, 1, 1));
end

% Other initializations for batch processing
stype = parse_type(ttype);  % Determine the size in bytes of the data type
dss = 2;                    % Spatial downsampling factor for visualization
dst = 2;                    % Temporal downsampling factor for visualization
nf = size(sigfn, 2);        % Number of frames in 'sigfn'
nsize = pixh * pixw * nf * stype * 6 / (dss ^ 2);  % Estimate the size of data to be processed
nbatch = batch_compute(nsize);                     % Compute the number of batches based on available memory
ebatch = ceil(nf / nbatch);                        % Number of frames per batch
idbatch = [1: ebatch: nf, nf + 1];                 % Indices to split data into batches
nbatch = length(idbatch) - 1;                      % Actual number of batches

% Define the filename for the TIFF stack
% Extract the path and base filename from 'fname'
[path_str, base_filename, ~] = fileparts(fname);

% Define the filename for the TIFF stack, including the full path
tiff_filename = fullfile(path_str, [base_filename, '_processed_data.tif']);

% Delete the existing file if it exists
if exist(tiff_filename, 'file')
    delete(tiff_filename);
end

% Initialize a frame counter for the TIFF stack
frame_counter = 0;

% Loop over batches to process and save frames
for ii = 1: nbatch
    % Load data based on motion correction flag
    if ismc
        % Use motion-corrected data from 'frame_all'
        dataraw = mraw.frame_all(1:dss:pixh, 1:dss:pixw, idbatch(ii):idbatch(ii+1)-1);
    else
        % Use preprocessed data without motion correction from 'frame_allt'
        dataraw = mraw.frame_allt(1:dss:pixh, 1:dss:pixw, idbatch(ii):idbatch(ii+1)-1);
    end
    % Load registered data (if applicable)
    datareg = mreg.reg(1:dss:pixh, 1:dss:pixw, idbatch(ii):idbatch(ii+1)-1);

    % Compute the processed data if possible
    if ~isempty(sigfn) && ~isempty(roifn)
        % Multiply spatial filters (roifn) by temporal signals (sigfn) to reconstruct the fluorescence signal
        datar = reshape(roifn * sigfn(:, idbatch(ii):idbatch(ii+1)-1), pixh, pixw, []);
        datar = datar(1:dss:end, 1:dss:end, :);  % Apply spatial downsampling
    else
        datar = [];  % No processed data available
    end

    % Loop over frames within the batch to save them to the TIFF stack
    num_frames_in_batch = size(dataraw, 3);
    for i = 1:dst:num_frames_in_batch
        frame_counter = frame_counter + 1;
        % Prepare the processed frame
        if ~isempty(datar)
            processed_frame = datar(:, :, i);
        else
            processed_frame = zeros(size(dataraw(:, :, i)));  % Use a blank frame if no data
        end

        % Convert the processed frame to uint16 format (assuming data ranges between 0 and 1)
        processed_frame_uint16 = uint16(processed_frame * 65535);

        % Write the processed frame to the TIFF stack
        if frame_counter == 1
            % For the first frame, overwrite any existing file
            imwrite(processed_frame_uint16, tiff_filename, 'WriteMode', 'overwrite', 'Compression', 'none');
        else
            % For subsequent frames, append to the TIFF stack
            imwrite(processed_frame_uint16, tiff_filename, 'WriteMode', 'append', 'Compression', 'none');
        end
    end
end

%% Create the Green-Hot Colormap %%
N = 256;  % Number of colors in the colormap
% Define color points: black -> green -> white
cmap_points = [
    0 0 0;   % Black
    0 1 0;   % Green
    1 1 1    % White
    ];
% Shift the midpoint towards lower intensities (e.g., 0.3) to enhance visibility of low-intensity values
color_positions = [0; 0.2; 1];  % Positions of the colors in the colormap (from 0 to 1)
% Interpolate to create the colormap
greenhot = interp1(color_positions, cmap_points, linspace(0, 1, N));

%%% Compute Max Value of 'datar' for Consistent Color Scaling %%%
max_datar = 0;

% Check if 'sigfn' and 'roifn' are not empty
if ~isempty(sigfn) && ~isempty(roifn)
    for ii = 1: nbatch
        % Compute 'datar_batch' for the current batch
        datar_batch = reshape(roifn * sigfn(:, idbatch(ii): idbatch(ii + 1) -1), pixh, pixw, []);
        datar_batch = datar_batch(1: dss: end, 1: dss: end, :);
        if ~isempty(datar_batch)
            max_datar = max(max_datar, max(datar_batch(:)));
        end
    end

    % Ensure 'max_datar' is valid
    if isempty(max_datar) || max_datar <= 0
        disp('No valid data for processed frames. Setting max_datar to 1.');
        max_datar = 1;  % Set to default positive value
    end
else
    disp('No cells found. Skipping processed data computation.');
    max_datar = 1;  % Set to default value
end

%% Make a Movie %%
% Define the filename for the movie, including the full path
fmovie = fullfile(path_str, [base_filename, '_demo_vid.avi']);

% Set up the video writer object
v = VideoWriter(fmovie);
v.FrameRate = Fsi_new;  % Set the frame rate to the new sampling rate
v.Quality = 100;        % Set the video quality to maximum
open(v)

%%% Generate the Movie %%%
figure(2)
set(gcf, 'Units', 'normalized', 'position', [0.2, 0.1, 0.6, 0.6])  % Adjust the figure size for better visibility

for ii = 1: nbatch
    % Load raw data based on motion correction flag
    if ismc
        % Use motion-corrected data from 'frame_all'
        dataraw = mraw.frame_all(1: dss: pixh, 1: dss: pixw, idbatch(ii): idbatch(ii + 1) - 1);
    else
        % Use preprocessed data without motion correction from 'frame_allt'
        dataraw = mraw.frame_allt(1: dss: pixh, 1: dss: pixw, idbatch(ii): idbatch(ii + 1) - 1);
    end
    % Load registered data (if applicable)
    datareg = mreg.reg(1: dss: pixh, 1: dss: pixw, idbatch(ii): idbatch(ii + 1) - 1);

    % Compute the processed data if possible
    if ~isempty(sigfn) && ~isempty(roifn)
        datar = reshape(roifn * sigfn(:, idbatch(ii): idbatch(ii + 1) - 1), pixh, pixw, []);
        datar = datar(1: dss: end, 1: dss: end, :);  % Apply spatial downsampling
    else
        datar = [];  % No processed data available
    end

    num_frames_in_batch = size(dataraw, 3);
    for i = 1: dst: num_frames_in_batch
        clf  % Clear the current figure

        %--- Plot Raw Data ---
        ax1 = subplot(1, 2, 1, 'align');
        imagesc(dataraw(:, :, i), [0, 1])
        axis off
        axis square
        title('Raw')
        colormap(ax1, 'gray')  % Use grayscale colormap for raw data

        %--- Plot Processed Data with Custom Colormap ---
        if ~isempty(datar) && max_datar > 0
            ax3 = subplot(1, 2, 2, 'align');
            imagesc(datar(:, :, i), [0, max_datar])  % Use consistent color scaling
            axis off
            axis square
            title('Processed')
            colormap(ax3, greenhot)  % Apply the custom colormap to this subplot
        else
            % No processed data to plot
            subplot(1, 2, 2, 'align');
            axis off
            axis square
            title('No Processed Data')
        end

        % Add a super title with the frame number
        sgtitle(['Frame #', num2str(i + idbatch(ii) - 1)])

        % Capture the current figure as a frame in the video
        movtmp = getframe(gcf);
        writeVideo(v, movtmp);
    end
end
close(v)  % Close the video file

%% Post-Processing (Optional) %%
if ifpost
    real_neuron_select  % Call a function for additional post-processing (if desired)
end


% %%% demo of the full MIN1PIPE %%%
% 
% clear all;
% 
% %% session-specific parameter initialization %%
% Fsi = 40;
% Fsi_new = 40; %%% no temporal downsampling %%%
% spatialr = 0.5; %%% no spatial downsampling %%%
% se = 9; %%% structure element for background removal %%%
% ismc = 0; %%% run movement correction %%%
% flag = 2; %%% use auto seeds selection; 2 if manual %%%
% %isvis = true; %%% do visualize %%%
% ifpost = false; %%% set true if want to see post-process %%%
% 
% %% main program %%
% [fname, frawname, fregname] = min1pipe(Fsi, Fsi_new, spatialr, se, ismc, flag); % old if you want to manually select parameters %
% %[fname, frawname, fregname] = min1pipe(Fsi, Fsi_new, [], [], ismc, flag);
% 
% %% plot some images %%
% load(fname)
% figure(1)
% clf
% %%% raw max %%%
% subplot(2, 3, 1, 'align')
% imagesc(imaxn)
% axis square
% title('Raw')
% 
% %%% neural enhanced before movement correction %%%
% subplot(2, 3, 2, 'align')
% imagesc(imaxy)
% axis square
% title('Before MC')
% 
% %%% neural enhanced after movement correction %%%
% subplot(2, 3, 3, 'align')
% imagesc(imax)
% axis square
% title('After MC')
% 
% %%% contour %%%
% subplot(2, 3, 4, 'align')
% plot_contour(roifn, sigfn, seedsfn, imax, pixh, pixw)
% axis square
% % if 0
% % %%% movement measurement %%%
% % subplot(2, 3, 5, 'align')
% % axis off
% % if ismc
% %     plot(raw_score); hold on; plot(corr_score); hold off;
% %     axis square
% %     title('MC Scores')
% % else
% %     title('MC skipped')
% % end
% % end
% 
% % Plotting traces
% subplot(2, 3, [5 6])
% if ~isempty(sigfn)
%     sigt = sigfn;
%     for i = 1: size(sigt, 1)
%         sigt(i, :) = normalize(sigt(i, :));
%     end
%     plot((sigt + (1: size(sigt, 1))')')
%     axis tight
%     title('Traces')
% else
%     % Display message or placeholder
%     axis off
%     title('No Traces to Display')
% end
% 
% 
% %% making a tiff stack
% load(fname)
% mraw = matfile(frawname);
% mreg = matfile(fregname);
% % Determine the data type
% if ismc
%     ttype = class(mraw.frame_all(1, 1, 1));
% else
%     ttype = class(mraw.frame_allt(1, 1, 1));
% end
% 
% % Other initializations
% stype = parse_type(ttype);
% dss = 2;
% dst = 2;
% nf = size(sigfn, 2);
% nsize = pixh * pixw * nf * stype * 6 / (dss ^ 2); %%% size of single %%%
% nbatch = batch_compute(nsize);
% ebatch = ceil(nf / nbatch);
% idbatch = [1: ebatch: nf, nf + 1];
% nbatch = length(idbatch) - 1;
% 
% % Define the filename for the TIFF stack
% % Extract the path and base filename from 'fname'
% [path_str, base_filename, ~] = fileparts(fname); % 'path_str' gives the folder, 'base_filename' gives the file name without extension
% 
% % Define the filename for the TIFF stack, including the full path
% tiff_filename = fullfile(path_str, [base_filename, '_processed_data.tif']);
% 
% 
% 
% 
% % Delete the existing file if it exists
% if exist(tiff_filename, 'file')
%     delete(tiff_filename);
% end
% 
% % Initialize a frame counter
% frame_counter = 0;
% 
% % Loop over batches
% for ii = 1: nbatch
% 
%     % Compute processed data if possible
%     if ~isempty(sigfn) && ~isempty(roifn)
%         datar = reshape(roifn * sigfn(:, idbatch(ii): idbatch(ii+1)-1), pixh, pixw, []);
%         datar = datar(1:dss:end, 1:dss:end, :);
%     else
%         datar = []; % No processed data
%     end
% 
% 
%     for i = 1: dst: size(dataraw, 3)
%         clf
%         %--- Plot Raw Data ---
%         ax1 = subplot(1, 2, 1, 'align');
%         imagesc(dataraw(:, :, i), [0, 1])
%         axis off
%         axis square
%         title('Raw')
%         colormap(ax1, 'gray')  % Use grayscale colormap for raw data
% 
%         %--- Plot Processed Data ---
%         if ~isempty(datar) && max_datar > 0
%             ax3 = subplot(1, 2, 2, 'align');
%             imagesc(datar(:, :, i), [0, max_datar])  % Use consistent color scaling
%             axis off
%             axis square
%             title('Processed')
%             colormap(ax3, greenhot)  % Apply the custom colormap to this subplot
%         else
%             % No processed data to plot
%             subplot(1, 2, 2, 'align');
%             axis off
%             axis square
%             title('No Processed Data')
%         end
% 
%         sgtitle(['Frame #', num2str(i + idbatch(ii) - 1)])
% 
%         % Capture the frame and write to video
%         movtmp = getframe(gcf);
%         writeVideo(v, movtmp);
%     end
% 
%     % Load data based on motion correction flag
%     if ismc
%         dataraw = mraw.frame_all(1:dss:pixh, 1:dss:pixw, idbatch(ii):idbatch(ii+1)-1);
%     else
%         dataraw = mraw.frame_allt(1:dss:pixh, 1:dss:pixw, idbatch(ii):idbatch(ii+1)-1);
%     end
%     datareg = mreg.reg(1:dss:pixh, 1:dss:pixw, idbatch(ii):idbatch(ii+1)-1);
% 
%     % Compute the processed data
%     datar = reshape(roifn * sigfn(:, idbatch(ii):idbatch(ii+1)-1), pixh, pixw, []);
%     datar = datar(1:dss:end, 1:dss:end, :);
% 
%     % Loop over frames within the batch
%     for i = 1:dst:size(datar, 3)
%         frame_counter = frame_counter + 1;
%         processed_frame = datar(:, :, i);
% 
%         % Convert the processed frame to uint16 format
%         % Assuming processed_frame is of type double and ranges between 0 and 1
%         processed_frame_uint16 = uint16(processed_frame * 65535);
% 
%         % Write the processed frame to the TIFF stack
%         if frame_counter == 1
%             % For the first frame, overwrite any existing file
%             imwrite(processed_frame_uint16, tiff_filename, 'WriteMode', 'overwrite', 'Compression', 'none');
%         else
%             % For subsequent frames, append to the TIFF stack
%             imwrite(processed_frame_uint16, tiff_filename, 'WriteMode', 'append', 'Compression', 'none');
%         end
%     end
% end
% 
% 
% 
% %% Create the green-hot colormap
% N = 256;  % Number of colors in the colormap
% % Define color points: black -> green -> white
% cmap_points = [
%     0 0 0;   % Black
%     0 1 0;   % Green
%     1 1 1    % White
%     ];
% color_positions = [0; 0.3; 1];  % Positions of the colors in the colormap (from 0 to 1)
% % Interpolate to create the colormap
% greenhot = interp1(color_positions, cmap_points, linspace(0, 1, N));
% 
% 
% %%% Compute batch %%%
% stype = parse_type(ttype);
% dss = 2;
% dst = 2;
% nf = size(sigfn, 2);
% nsize = pixh * pixw * nf * stype * 6 / (dss ^ 2); %%% Size of single %%%
% nbatch = batch_compute(nsize);
% ebatch = ceil(nf / nbatch);
% idbatch = [1: ebatch: nf, nf + 1];
% nbatch = length(idbatch) - 1;
% 
% %%% Compute max value of 'datar' for consistent color scaling %%%
% max_datar = 0;
% for ii = 1: nbatch
%     datar_batch = reshape(roifn * sigfn(:, idbatch(ii): idbatch(ii + 1) -1), pixh, pixw, []);
%     datar_batch = datar_batch(1: dss: end, 1: dss: end, :);
%     max_datar = max(max_datar, max(datar_batch(:)));
% end
% 
% 
% 
% %% make a movie %%
% load(fname)
% % Extract the path and base filename from 'fname'
% [path_str, base_filename, ~] = fileparts(fname); % 'path_str' gives the folder, 'base_filename' gives the file name without extension
% 
% % Define the filename for the movie, including the full path
% fmovie = fullfile(path_str, [base_filename, '_demo_vid.avi']);
% 
% v = VideoWriter(fmovie);
% v.FrameRate = Fsi_new;
% v.Quality = 100;
% open(v)
% 
% %%% compute batch %%%
% 
% % Other initializations
% stype = parse_type(ttype);
% dss = 2;
% dst = 2;
% nf = size(sigfn, 2);
% nsize = pixh * pixw * nf * stype * 6 / (dss ^ 2); %%% size of single %%%
% nbatch = batch_compute(nsize);
% ebatch = ceil(nf / nbatch);
% idbatch = [1: ebatch: nf, nf + 1];
% nbatch = length(idbatch) - 1;
% 
% %%% make movie %%%
% 
% %%% Make movie %%%
% figure(2)
% set(gcf, 'Units', 'normalized', 'position', [0.2, 0.1, 0.6, 0.6])  % Adjusted for better visibility
% for ii = 1: nbatch
%     if ismc
%         dataraw = mraw.frame_all(1: dss: pixh, 1: dss: pixw, idbatch(ii): idbatch(ii + 1) - 1);
%     else
%         dataraw = mraw.frame_allt(1: dss: pixh, 1: dss: pixw, idbatch(ii): idbatch(ii + 1) - 1);
%     end
%     datareg = mreg.reg(1: dss: pixh, 1: dss: pixw, idbatch(ii): idbatch(ii + 1) - 1);
%     datar = reshape(roifn * sigfn(:, idbatch(ii): idbatch(ii + 1) - 1), pixh, pixw, []);
%     datar = datar(1: dss: end, 1: dss: end, :);
%     for i = 1: dst: size(dataraw, 3)
%         clf
%         %--- Plot Raw Data ---
%         ax1 = subplot(1, 2, 1, 'align');
%         imagesc(dataraw(:, :, i), [0, 1])
%         axis off
%         axis square
%         title('Raw')
%         colormap(ax1, 'gray')  % Use grayscale colormap for raw data
% 
%         % %--- Plot After Motion Correction ---
%         % ax2 = subplot(1, 3, 2, 'align');
%         % imagesc(datareg(:, :, i), [0, 1])
%         % axis off
%         % axis square
%         % title('After MC')
%         % colormap(ax2, 'gray')  % Use grayscale colormap for registered data
% 
%         %--- Plot Processed Data with Custom Colormap ---
%         ax3 = subplot(1, 2, 2, 'align');
%         imagesc(datar(:, :, i), [0, max_datar])  % Use consistent color scaling
%         axis off
%         axis square
%         title('Processed')
%         colormap(ax3, greenhot)  % Apply the custom colormap to this subplot
% 
%         sgtitle(['Frame #', num2str(i + idbatch(ii) - 1)])
% 
%         % Capture the frame and write to video
%         movtmp = getframe(gcf);
%         writeVideo(v, movtmp);
%     end
% end
% close(v)
% 
% %% post-process %%
% if ifpost
%     real_neuron_select
% end