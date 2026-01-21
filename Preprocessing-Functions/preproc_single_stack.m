%% 1. Load files into single 4d array
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

%% 2. Masking muscle of interest

anat_image = dti_all_unreg(:, :, :, 1);

% Define the slice range (process all slices)
% Input format is [Start_Slice, End_Slice]
slice_range = [1, slices]; 

% Define options (leave empty to skip complex features)
alt_mask_size = []; % We don't need to resize since we are working in DTI space
fv_options = [];    % Skip advanced visualization for now

fprintf("Opening Manual Segmentation Tool...\n");
fprintf("Instructions:\n");
fprintf("  1. Left-click to draw points around the muscle.\n");
fprintf("  2. Right-click inside the shape and select 'Create Mask' to finish a slice.\n");
fprintf("  3. Repeat for all slices.\n");

% Call the toolbox function
[pd_mask, ~] = define_muscle(anat_image, slice_range, alt_mask_size, fv_options);