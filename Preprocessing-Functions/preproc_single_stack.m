%% 1.1 Load diffusion files into single 4d array
series_path = uigetdir(pwd, 'Select the Folder containing the DTI DICOMs series');

cd(series_path);
file_list = dir(fullfile(series_path, '*.dcm'));
num_files = length(file_list);

[~, idx] = sort(str2double(regexp({file_list.name}, '\d+', 'match', 'once'))); %sort file list
file_list = file_list(idx);

% --- EXTRACT DTI GEOMETRY FROM FIRST DICOM HEADER ---
temp_info = dicominfo(file_list(1).name);

if isfield(temp_info, 'SliceThickness')
    dti_slice_thick = temp_info.SliceThickness;
elseif isfield(temp_info, 'SpacingBetweenSlices')
    dti_slice_thick = temp_info.SpacingBetweenSlices;
else
    dti_slice_thick = 7; % Fallback
    fprintf('WARNING: Could not find DTI Slice Thickness. Defaulting to 7mm.\n');
end

if isfield(temp_info, 'PixelSpacing')
    dti_fov = temp_info.PixelSpacing(1) * double(temp_info.Rows); 
else
    dti_fov = 192; % Fallback
    fprintf('WARNING: Could not find DTI Pixel Spacing. Defaulting FOV to 192mm.\n');
end

dti_dims = [dti_fov, dti_slice_thick];
fprintf('DTI Geometry calculated: FOV = %.1f mm, Slice Thickness = %.1f mm\n', dti_fov, dti_slice_thick);
% ----------------------------------------------------

temp_vol = squeeze(double(dicomread(file_list(1).name))); %get dimensions from first file
[rows, cols, slices] = size(temp_vol);
dti_all_unreg = zeros(rows, cols, slices, num_files);

bvect = zeros(num_files, 3); % Create containers for Gradient Directions and B-Values
bval_list = zeros(num_files, 1);

for i = 1:num_files
    currentFile = file_list(i).name;

    vol_data = squeeze(double(dicomread(currentFile)));
    dti_all_unreg(:, :, :, i) = vol_data; % Store the volume data in the 4D array

    info = dicominfo(currentFile);
    seq = info.PerFrameFunctionalGroupsSequence.Item_1.MRDiffusionSequence.Item_1;
    bval_list(i) = seq.DiffusionBValue;
    if isfield(seq, 'DiffusionGradientDirectionSequence')
        vec = seq.DiffusionGradientDirectionSequence.Item_1.DiffusionGradientOrientation;
        bvect(i,:) = vec;
    else % If no direction exists, it's likely the b=0 image
        bvect(i,:) = [0 0 0];
        fprintf("Found b=0 image!\n");
    end
    fprintf('Adding file %d: %s\n', i, currentFile);
end
fprintf('Averaging b=0 images...\n');

b0_inds = find(bval_list <= 10); %isolate b=0 maps from valid directions
dwi_inds = find(bval_list > 10);
fprintf("Found %d b=0 images and %d diffusion directions.\n", length(b0_inds), length(dwi_inds));

b0_avg = mean(dti_all_unreg(:, :, :, b0_inds), 4); %average along the 4th dimension (time/direction dimension)

dwi_data = dti_all_unreg(:, :, :, dwi_inds); % reconstruct the final 4D array
dti_all_unreg = cat(4, b0_avg, dwi_data); % structure: [Average_b0, DWI_1, DWI_2, ... DWI_N]
bvect = bvect(dwi_inds, :);

% Extract the single diffusion b-value for tensor calculation later
bval = max(bval_list); 

% Default to false (will be overwritten if Section 1.2 is run successfully)
using_anat = false;

fprintf("Done with loading!\n");

% Clean up workspace, explicitly keeping our new variables (dti_dims, bval, and using_anat)
clear series_path file_list num_files idx temp_info dti_slice_thick dti_fov temp_vol i currentFile vol_data info seq vec bval_list b0_inds dwi_inds b0_avg dwi_data

% Save the DTI volume, dimensions, gradient vectors, the b-value, and the anatomical flag
save('dti_all_unreg.mat', 'dti_all_unreg', 'dti_dims', 'using_anat', '-v7.3');
save('bvect.mat', 'bvect', 'bval');

%% 1.2 Load anatomical volume (if available)
[anat_file, anat_dir] = uigetfile('*.dcm', 'Select the anatomical volume (if available, otherwise click cancel)');

if isequal(anat_file, 0)
    fprintf("No anatomical volume selected.\n");
else
    fprintf("Loading anatomical volume...\n");
    full_anat_path = fullfile(anat_dir, anat_file);
    anat_vol = squeeze(double(dicomread(full_anat_path)));
    
    % --- EXTRACT ANATOMICAL GEOMETRY FROM DICOM HEADER ---
    anat_info = dicominfo(full_anat_path);
    
    if isfield(anat_info, 'SliceThickness')
        anat_slice_thick = anat_info.SliceThickness;
    else
        anat_slice_thick = 7; % Fallback if missing
        fprintf('WARNING: Could not find Anatomical Slice Thickness. Defaulting to 7mm.\n');
    end
    
    if isfield(anat_info, 'PixelSpacing')
        anat_fov = anat_info.PixelSpacing(1) * double(anat_info.Rows);
    else
        anat_fov = 192; % Fallback if missing
        fprintf('WARNING: Could not find Anatomical Pixel Spacing. Defaulting FOV to 192mm.\n');
    end
    
    anat_dims = [anat_fov, anat_slice_thick];
    fprintf('Anatomical Geometry calculated: FOV = %.1f mm, Slice Thickness = %.1f mm\n', anat_fov, anat_slice_thick);
    % -----------------------------------------------------

    % Resize anatomical volume to match DTI dimensions
    fprintf('Resizing Anatomical data to match DTI geometry...\n');
    [target_rows, target_cols, target_slices, ~] = size(dti_all_unreg);
    
    % Use 'linear' (trilinear) interpolation to prevent Gibbs ringing artifacts
    anat_vol = imresize3(anat_vol, [target_rows, target_cols, target_slices], 'linear'); 
    
    % CRITICAL: Once resized, the Anatomical volume now physically occupies the DTI spatial matrix!
    if exist('dti_dims', 'var')
        anat_dims = dti_dims; 
    else
        fprintf('WARNING: dti_dims variable not found in workspace. Make sure Section 1.1 extracted it!\n');
    end
    
    % Successfully loaded and resized an anatomical image, so overwrite the flag
    using_anat = true;
    
    fprintf('Resizing complete!\n');
end

clear anat_file anat_dir target_rows target_cols target_slices full_anat_path anat_info anat_slice_thick anat_fov

%% 2.1 Masking muscle of interest manually with volume segmenter

% Use a temporary reference volume so we don't trick Section 3
if exist('using_anat', 'var') && using_anat && exist('anat_vol', 'var')
    fprintf("Using existing anatomical volume for reference.\n");
    ref_vol = anat_vol;
else
    fprintf("No anatomical volume found. Using average b=0 image for anatomical reference.\n");
    ref_vol = dti_all_unreg(:, :, :, 1); 
end

% Open the app with this data loaded
volumeSegmenter(ref_vol)

clear ref_vol % Keep workspace clean

%% 2.2 Otsu's Method for automatic masking

fprintf("Automatically generating mask...\n");

% Use a temporary reference volume
if exist('using_anat', 'var') && using_anat && exist('anat_vol', 'var')
    fprintf("Using existing anatomical volume for reference.\n");
    ref_vol = anat_vol;
else
    fprintf("No anatomical volume found. Using average b=0 image for anatomical reference.\n");
    ref_vol = dti_all_unreg(:, :, :, 1); 
end

% Explicitly define dimensions so the loop doesn't crash
[rows, cols, slices] = size(ref_vol);
pd_mask = zeros(rows, cols, slices);

for z = 1:slices
    loop_image = ref_vol(:,:,z);
    
    % Normalize the slice
    if max(loop_image(:)) > 0
        loop_image = loop_image / max(loop_image(:)); 
    end
    
    % Form mask using Otsu's threshold
    loop_mask = zeros(rows, cols);
    threshold = graythresh(loop_image); 
    loop_mask(loop_image > threshold) = 1;
    
    % Morphological cleanup (fill holes, erode edges to remove skin/subcutaneous fat)
    % fills small gap but not large (like bone)
    loop_mask = bwareaopen(logical(loop_mask), 500);
    inverted_mask = ~loop_mask;
    inverted_mask = bwareaopen(inverted_mask, 200); 
    loop_mask = ~inverted_mask;

    loop_mask = bwmorph(loop_mask, 'erode', 2); % You may need to adjust the '2' based on the dataset
    
    pd_mask(:,:,z) = loop_mask;
end

clear z loop_image threshold loop_mask ref_vol rows cols slices inverted_mask
save('pd_mask.mat', 'pd_mask', '-v7.3');
fprintf("Mask succesfully generated and saved!\n");

%% 2.3 Upload mask from mat file
pd_mask = load(uigetfile('*.mat', 'Select the mask'));

% Ensure data is loaded in as a matrix, not a struct
if isstruct(pd_mask) 
    fields = fieldnames(pd_mask);
    pd_mask = pd_mask.(fields{1});
end

fprintf("Mask succesfully loaded!\n");

%% 2.4 Upload mask from nifti file (if this fails try loading from mat file)

%TODO: add nifti loading functionality here (using niftiread and niftiinfo)

%% 2.5 define_muscle built in masking method

fprintf("Starting manual segmentation using define_muscle...\n");

% Use a temporary reference volume
if exist('using_anat', 'var') && using_anat && exist('anat_vol', 'var')
    fprintf("Using existing anatomical volume for reference.\n");
    ref_vol = anat_vol;
else
    fprintf("No anatomical volume found. Using average b=0 image for anatomical reference.\n");
    ref_vol = dti_all_unreg(:, :, :, 1); 
end

% determine total slices available
total_slices = size(ref_vol, 3);

% prompt user for the slice range they want to segment
prompt = {sprintf('Enter starting slice (1-%d):', total_slices), ...
          sprintf('Enter ending slice (1-%d):', total_slices)};
dlgtitle = 'Slice Range for define_muscle';
dims = [1 50];
definput = {'1', num2str(total_slices)};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    fprintf("Segmentation canceled by user.\n");
else
    slices_to_segment = [str2double(answer{1}), str2double(answer{2})];

    % set up fiber_visualizer options dynamically
    if exist('anat_dims', 'var')
        current_dims = anat_dims;
    else
        current_dims = dti_dims;
    end

    fv_options.anat_dims = current_dims; 
    fv_options.anat_slices = slices_to_segment(1):2:slices_to_segment(2); 
    
    fv_options.plot_mesh = 0;   
    fv_options.plot_mask = 1;   
    fv_options.plot_fibers = 0; 
    
    fv_options.mask_size = [size(ref_vol, 1) size(ref_vol, 2)];
    fv_options.mask_dims = current_dims; % Automatically sets to [FOV SliceThickness]
    fv_options.mask_color = [1 0 0];

    fprintf("Interactive window opening. Follow the define_muscle instructions.\n");
    
    % call the function using our temporary reference volume
    [pd_mask, ~] = define_muscle(ref_vol, slices_to_segment, [], fv_options);
    
    save('pd_mask.mat', 'pd_mask', '-v7.3');
    fprintf("define_muscle complete! Variable 'pd_mask' saved to disk.\n");
end

clear total_slices prompt dlgtitle dims definput answer slices_to_segment fv_options ref_vol current_dims

%% 3.1 Demons registration method

fprintf('Starting Registration (this may take time)....\n');

% 1. Ensure the mask is a clean matrix
if isstruct(pd_mask) 
    fields = fieldnames(pd_mask);
    pd_mask = pd_mask.(fields{1});
    fprintf("Fixed pd_mask struct. Now it is a matrix.\n");
end

% setup volumes and pre-allocate
[rows, cols, slcs, total_vols] = size(dti_all_unreg);
dti_all_reg = zeros(size(dti_all_unreg));

% choose ground truth reference image
if exist('using_anat', 'var') && using_anat && exist('anat_vol', 'var')
    fprintf('Using Anatomical volume as fixed reference for registration.\n');
    fixed_vol = anat_vol;
else
    fprintf('Using average b=0 volume as fixed reference for registration.\n');
    fixed_vol = dti_all_unreg(:, :, :, 1);
    using_anat = false; % enforce false just in case
    dti_all_reg(:, :, :, 1) = dti_all_unreg(:, :, :, 1); % b=0 is perfectly aligned with itself
end

% registration loop slice by slice
fprintf('Starting Slice-by-Slice Demons Registration...\n');

for v = 1:total_vols 
    
    % Skip registering the b=0 image to itself if we aren't using an anatomical reference
    if ~using_anat && v == 1
        continue; 
    end

    fprintf('  Registering Volume %d of %d...\n', v, total_vols);
    
    for z = 1:slcs 
        % Grab the slices and current mask
        fixed_slice = fixed_vol(:, :, z); 
        moving_slice = dti_all_unreg(:, :, z, v);
        current_mask = pd_mask(:, :, z);

        % Apply mask BEFORE normalization to isolate the muscle contrast
        fixed_masked = fixed_slice .* current_mask;
        moving_masked = moving_slice .* current_mask;

        % Safely normalize (Skip registration for this slice if the mask is empty)
        fixed_max = max(fixed_masked(:));
        moving_max = max(moving_masked(:));
        
        if fixed_max > 0 && moving_max > 0
            fixed_norm = fixed_masked / fixed_max; 
            moving_norm = moving_masked / moving_max;
        else
            dti_all_reg(:, :, z, v) = moving_slice; % Copy original and skip
            continue;
        end

        % THE CONTRAST HACK: Apply only when registering b=0 to Anatomical
        if using_anat && v == 1 
            moving_norm(moving_norm > 0.5) = 0.5; % Flatten out bright blood vessels
            moving_norm = moving_norm * 2;        % Scale up to match anatomical brightness
        end

        % Calculate displacement field and apply the warp to the UNMASKED original slice
        [disp_field, ~] = imregdemons(moving_norm, fixed_norm, [500 400 200], 'AccumulatedFieldSmoothing', 1.3, 'DisplayWaitbar', false);
        dti_all_reg(:, :, z, v) = imwarp(moving_slice, disp_field); 
    end
end

save('dti_registered.mat', 'dti_all_reg', 'pd_mask', 'bvect', '-v7.3');
fprintf('Registration Complete! Variable "dti_all_reg" saved.\n');

clear fields rows cols slcs total_vols fixed_vol v z fixed_slice moving_slice current_mask fixed_masked moving_masked fixed_max moving_max fixed_norm moving_norm disp_field

%% 3.2 FSL Eddy based registraton method

%TODO: base of off exhisting implementation in eddy_correct.m

%% 4.1 Denoising with anisotropic smoothing

%TODO: Implement based of script provided

%% 4.2 Threshold based principle component analysis (PCA) denoising

%TODO: Implement based of script provided