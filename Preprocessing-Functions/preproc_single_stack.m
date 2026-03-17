%% 1.1 Load diffusion files into single 4d array
series_path = uigetdir(pwd, 'Select the Folder containing the DTI DICOMs series');

cd(series_path);
file_list = dir(fullfile(series_path, '*.dcm'));
num_files = length(file_list);

[~, idx] = sort(str2double(regexp({file_list.name}, '\d+', 'match', 'once'))); %sort file list
file_list = file_list(idx);

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

fprintf("Done with loading!\n");

clear series_path file_list num_files idx temp_vol i currentFile vol_data info seq vec bval_list b0_inds dwi_inds b0_avg dwi_data

%% 1.2 Load anatomical volume (if available)
[anat_file, anat_dir] = uigetfile('*.dcm', 'Select the anatomical volume (if available, otherwise click cancel)');

if isequal(anat_file, 0)
    fprintf("No anatomical volume selected.\n");
else
    fprintf("Loading anatomical volume...\n");
    anat_vol = squeeze(double(dicomread(fullfile(anat_dir, anat_file))));
    
    %resize anatomical volume to match DTI dimensions (if necessary)
    fprintf('Resizing Anatomical data to match DTI geometry...\n')

    [target_rows, target_cols, target_slices, ~] = size(dti_all_unreg)
    anat_vol = imresize3(anat_vol, [target_rows, target_cols, target_slices]); %weighted avg of 8 closest pixels in 3d space to resize the anatomical volume to match the DTI dimensions

    printf('Resizing complete!\n');
    clear anat_filename anat_filepath full_anat_path anat_orig target_rows target_cols target_slices
end

%% 2.1 Masking muscle of interest manually

% Check for existing anatomical volume, otherwise use average b=0 image
if ~exist('anat_vol', 'var')
    fprintf("No anatomical volume found. Using average b=0 image for anatomical reference.\n");
    anat_vol = dti_all_unreg(:, :, :, 1); 
else
    fprintf("Using existing anatomical volume for reference.\n")
end
% Open the app with this data loaded
volumeSegmenter(anat_vol)

%% 2.2 Otsu's Method for automatic masking

fprintf("Automatically generating mask...\n");

% Check for existing anatomical volume, otherwise use average b=0 image
if ~exist('anat_vol', 'var')
    fprintf("No anatomical volume found. Using average b=0 image for anatomical reference.\n");
    anat_vol = dti_all_unreg(:, :, :, 1); 
else
    fprintf("Using existing anatomical volume for reference.\n");
end

pd_mask = zeros(size(anat_vol));

for z = 1:slices
    loop_image = anat_vol(:,:,z);
    
    % Normalize the slice
    if max(loop_image(:)) > 0
        loop_image = loop_image / max(loop_image(:)); 
    end
    
    % Form mask using Otsu's threshold
    loop_mask = zeros(rows, cols);
    threshold = graythresh(loop_image); 
    loop_mask(loop_image > threshold) = 1;
    
    % Morphological cleanup (fill holes, erode edges to remove skin/subcutaneous fat)
    loop_mask = imfill(loop_mask, 'holes');
    loop_mask = bwmorph(loop_mask, 'erode', 2); % You may need to adjust the '2' based on the dataset
    
    pd_mask(:,:,z) = loop_mask;
end

%% 2.3 Upload mask from mat file
pd_mask = load(uigetfile('Select the mask'));

fprintf("Mask succesfully loaded!\n");

%% 2.4 Upload mask from nifti file (if this fails try loading from mat file)\

%TODO: add nifti loading functionality here (using niftiread and niftiinfo)

%% 3.1 Demons registration method

fprintf('Starting Registration (this may take time)...\n');

if isstruct(pd_mask) % ensure data loaded in as matrix
    fields = fieldnames(pd_mask);
    pd_mask = pd_mask.(fields{1});
    fprintf("Fixed pd_mask struct. Now it is a matrix.\n");
end

fprintf('Starting Registration (Slice-by-Slice using Demons)...\n');

dti_all_reg = zeros(size(dti_all_unreg));
[rows, cols, slcs, total_vols] = size(dti_all_unreg);

dti_all_reg(:, :, :, 1) = dti_all_unreg(:, :, :, 1); % volume 1 is our fixed reference b=0

for v = 2:total_vols %for each volume
    fprintf('  Registering Volume %d of %d...\n', v, total_vols);
    for z = 1:slcs % for each slice in each volume
        fixed_slice = dti_all_unreg(:, :, z, 1); % get reference point and slice that is to be corrected
        moving_slice = dti_all_unreg(:, :, z, v);
        fixed_norm = fixed_slice / max(fixed_slice(:)); % normalize
        moving_norm = moving_slice / max(moving_slice(:));

        % calculate displacment field
        [disp_field, ~] = imregdemons(moving_norm, fixed_norm, [500 400 200], 'AccumulatedFieldSmoothing', 1.3, 'DisplayWaitbar', false);

        dti_all_reg(:, :, z, v) = imwarp(moving_slice, disp_field); % apply the warp to the slice
    end
end

save('dti_registered.mat', 'dti_all_reg', 'pd_mask', 'bvect', '-v7.3');

fprintf('Registration Complete! Variable "dti_all_reg" created.\n');

%% 3.2 FSL Eddy based registraton method

%TODO: base of off exhisting implementation in eddy_correct.m