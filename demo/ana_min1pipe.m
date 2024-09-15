function ana_min1pipe()
    % Main function for processing, visualizing, and saving results.
    
    %% Session-specific parameter initialization %%
    Fsi = 40;
    Fsi_new = 40; % No temporal downsampling
    spatialr = 1; % No spatial downsampling
    se = 9; % Structure element for background removal
    ismc = 0; % Run movement correction
    flag = 2; % Use auto seeds selection; 2 if manual

    %% Run the main pipeline %%
    [file_name_to_save, filename_raw, filename_reg] = min1pipe(Fsi, Fsi_new, spatialr, se, ismc, flag);
    
    %% Plot results %%
    plot_results(fname, imaxn, imaxy, imax, sigfn, roifn, seedsfn, pixh, pixw, ismc);

    %% Save processed data to a TIFF stack %%
    save_tiff_stack(fname, frawname, fregname, sigfn, roifn, ismc, pixh, pixw);

    %% Create video from the processed data %%
    create_video(fname, frawname, fregname, sigfn, roifn, ismc, pixh, pixw, Fsi_new);

    %% Post-processing (if needed) %%
    % If post-processing is enabled, call the post-processing function.
    ifpost = true; % Set this flag if post-processing is needed
    if ifpost
        real_neuron_select(); % Post-process neuron selection
    end
end

%% Plot Results %%
function plot_results(fname, imaxn, imaxy, imax, sigfn, roifn, seedsfn, pixh, pixw, ismc)
    figure(1)
    clf
    
    % Plot raw max image
    subplot(2, 3, 1, 'align')
    imagesc(imaxn)
    axis square
    title('Raw')

    % Plot neural enhanced before movement correction
    subplot(2, 3, 2, 'align')
    imagesc(imaxy)
    axis square
    title('Before MC')

    % Plot neural enhanced after movement correction
    subplot(2, 3, 3, 'align')
    imagesc(imax)
    axis square
    title('After MC')

    % Plot contours
    subplot(2, 3, 4, 'align')
    plot_contour(roifn, sigfn, seedsfn, imax, pixh, pixw)
    axis square

    % Plot identified traces
    subplot(2, 3, [5 6])
    sigt = normalize_traces(sigfn);
    plot((sigt + (1: size(sigt, 1))')')
    axis tight
    title('Traces')
end

%% Normalize Traces %%
function sigt = normalize_traces(sigfn)
    sigt = sigfn;
    for i = 1:size(sigfn, 1)
        sigt(i, :) = normalize(sigfn(i, :));
    end
end

%% Save Processed Data as TIFF Stack %%
function save_tiff_stack(fname, frawname, fregname, sigfn, roifn, ismc, pixh, pixw)
    % Load necessary data
    mraw = matfile(frawname);
    mreg = matfile(fregname);

    % Determine data type and prepare batch processing
    ttype = determine_data_type(mraw, ismc);
    [idbatch, nbatch] = setup_batches(ttype, pixh, pixw, sigfn);

    % Define the filename for the TIFF stack
    tiff_filename = [fname(1: find(fname == filesep, 1, 'last')), 'processed_data.tif'];
    if exist(tiff_filename, 'file'), delete(tiff_filename); end

    % Save each frame to the TIFF stack
    frame_counter = 0;
    for ii = 1:nbatch
        [dataraw, datareg] = load_batch(mraw, mreg, idbatch, ii, ismc, pixh, pixw);
        datar = process_data(roifn, sigfn, idbatch, ii, pixh, pixw);
        frame_counter = write_tiff_frames(datar, tiff_filename, frame_counter);
    end
end

%% Create Video %%
function create_video(fname, frawname, fregname, sigfn, roifn, ismc, pixh, pixw, Fsi_new)
    % Load necessary data
    mraw = matfile(frawname);
    mreg = matfile(fregname);

    % Determine data type and prepare batch processing
    ttype = determine_data_type(mraw, ismc);
    [idbatch, nbatch] = setup_batches(ttype, pixh, pixw, sigfn);

    % Prepare video file
    fmovie = [fname(1: find(fname == filesep, 1, 'last')), 'demo_vid.avi'];
    v = VideoWriter(fmovie);
    v.FrameRate = Fsi_new;
    v.Quality = 100;
    open(v)

    % Create a colormap
    greenhot = create_greenhot_colormap();

    % Compute max value for consistent color scaling
    max_datar = compute_max_datar(roifn, sigfn, idbatch, pixh, pixw);

    % Generate the video
    for ii = 1:nbatch
        [dataraw, datareg] = load_batch(mraw, mreg, idbatch, ii, ismc, pixh, pixw);
        datar = process_data(roifn, sigfn, idbatch, ii, pixh, pixw);
        save_video_frames(datar, dataraw, greenhot, v, max_datar, ii);
    end
    close(v)
end

%% Helper Functions for Video and TIFF Generation %%
function ttype = determine_data_type(mraw, ismc)
    if ismc
        ttype = class(mraw.frame_all(1, 1, 1));
    else
        ttype = class(mraw.frame_allt(1, 1, 1));
    end
end

function [idbatch, nbatch] = setup_batches(ttype, pixh, pixw, sigfn)
    stype = parse_type(ttype);
    dss = 2;
    dst = 2;
    nf = size(sigfn, 2);
    nsize = pixh * pixw * nf * stype * 6 / (dss ^ 2); % Size of single
    nbatch = batch_compute(nsize);
    idbatch = [1: ceil(nf / nbatch): nf, nf + 1];
end

function [dataraw, datareg] = load_batch(mraw, mreg, idbatch, ii, ismc, pixh, pixw)
    if ismc
        dataraw = mraw.frame_all(1:2:pixh, 1:2:pixw, idbatch(ii):idbatch(ii+1)-1);
    else
        dataraw = mraw.frame_allt(1:2:pixh, 1:2:pixw, idbatch(ii):idbatch(ii+1)-1);
    end
    datareg = mreg.reg(1:2:pixh, 1:2:pixw, idbatch(ii):idbatch(ii+1)-1);
end

function datar = process_data(roifn, sigfn, idbatch, ii, pixh, pixw)
    datar = reshape(roifn * sigfn(:, idbatch(ii):idbatch(ii+1)-1), pixh, pixw, []);
    datar = datar(1:2:end, 1:2:end, :);
end

function frame_counter = write_tiff_frames(datar, tiff_filename, frame_counter)
    for i = 1:size(datar, 3)
        frame_counter = frame_counter + 1;
        processed_frame_uint16 = uint16(datar(:, :, i) * 65535); % Convert to uint16
        if frame_counter == 1
            imwrite(processed_frame_uint16, tiff_filename, 'WriteMode', 'overwrite', 'Compression', 'none');
        else
            imwrite(processed_frame_uint16, tiff_filename, 'WriteMode', 'append', 'Compression', 'none');
        end
    end
end

function save_video_frames(datar, dataraw, greenhot, v, max_datar, ii)
    for i = 1:size(datar, 3)
        clf
        ax1 = subplot(1, 2, 1, 'align');
        imagesc(dataraw(:, :, i), [0, 1])
        axis off
        axis square
        colormap(ax1, 'gray')

        ax2 = subplot(1, 2, 2, 'align');
        imagesc(datar(:, :, i), [0, max_datar])
        axis off
        axis square
        colormap(ax2, greenhot)

        sgtitle(['Frame #', num2str(i + ii - 1)])
        movtmp = getframe(gcf);
        writeVideo(v, movtmp);
    end
end

function greenhot = create_greenhot_colormap()
    N = 256;
    cmap_points = [0 0 0; 0 1 0; 1 1 1];
    color_positions = [0; 0.5; 1];
    greenhot = interp1(color_positions, cmap_points, linspace(0, 1, N));
end

function max_datar = compute_max_datar(roifn, sigfn, idbatch, pixh, pixw)
    max_datar = 0;
    for ii = 1:length(idbatch)-1
        datar_batch = reshape(roifn * sigfn(:, idbatch(ii):idbatch(ii+1)-1), pixh, pixw, []);
        datar_batch = datar_batch(1:2:end, 1:2:end, :);
        max_datar = max(max_datar, max(datar_batch(:)));
    end
end
