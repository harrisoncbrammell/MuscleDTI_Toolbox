%% 1.1 Load diffusion files into single 4d array
series_path = uigetdir(pwd, 'Select the Folder containing the DTI DICOMs series');
output_dir  = uigetdir(pwd, 'Select output folder to save processed files');

cd(series_path);
file_list = dir(fullfile(series_path, '*.dcm'));
num_files = length(file_list);

[~, idx] = sort(str2double(regexp({file_list.name}, '\d+', 'match', 'once'))); %sort file list
file_list = file_list(idx);

% --- EXTRACT DTI GEOMETRY FROM FIRST DICOM HEADER ---
% Siemens Enhanced MR DICOMs (one file per volume) store geometry per-frame,
% not at the top level. PixelSpacing and SliceThickness live inside
% PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.
temp_info = dicominfo(file_list(1).name);

try
    per_frame = temp_info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1;
    dti_pixel_spacing = double(per_frame.PixelSpacing(1)); % in-plane voxel size (mm)
    % Use SpacingBetweenSlices (center-to-center) rather than SliceThickness
    % (slab thickness) for spatial filtering — this is the true inter-voxel distance.
    if isfield(per_frame, 'SpacingBetweenSlices')
        dti_slice_spacing = double(per_frame.SpacingBetweenSlices);
    else
        dti_slice_spacing = double(per_frame.SliceThickness);
        fprintf('WARNING: SpacingBetweenSlices missing, using SliceThickness = %.1f mm.\n', dti_slice_spacing);
    end
catch
    dti_pixel_spacing = 1.5; % Fallback
    dti_slice_spacing = 7.0; % Fallback
    fprintf('WARNING: Could not read per-frame geometry. Defaulting to %.1f mm in-plane, %.1f mm slice spacing.\n', dti_pixel_spacing, dti_slice_spacing);
end

% dti_dims = [in-plane pixel size, slice spacing] — used as voxel resolution
% by aniso4D_smoothing. NOT the FOV.
dti_dims = [dti_pixel_spacing, dti_slice_spacing];
fprintf('DTI Geometry: in-plane = %.4f mm, slice spacing = %.1f mm, matrix = %d x %d\n', ...
    dti_pixel_spacing, dti_slice_spacing, temp_info.Rows, temp_info.Columns);
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

b0_inds  = find(bval_list <= 10);
dwi_inds = find(bval_list > 10);
fprintf("Found %d b=0 images and %d diffusion directions.\n", length(b0_inds), length(dwi_inds));

% Report all unique b-values — multi-shell acquisitions require choosing one
% shell for tensor estimation (signal2tensor2 takes a single b-value).
unique_bvals = unique(bval_list(dwi_inds));
fprintf('Unique DWI b-values found: ');
fprintf('%.0f ', unique_bvals);
fprintf('\n');
if length(unique_bvals) > 1
    fprintf('WARNING: Multi-shell acquisition detected (%d shells).\n', length(unique_bvals));
    fprintf('signal2tensor2 requires a single b-value. Only the highest shell (b=%.0f) will be\n', max(unique_bvals));
    fprintf('used for tensor estimation. All shells are kept in dti_all_unreg for registration.\n');
end

b0_avg   = mean(dti_all_unreg(:, :, :, b0_inds), 4);
dwi_data = dti_all_unreg(:, :, :, dwi_inds);
dti_all_unreg = cat(4, b0_avg, dwi_data); % [b0_avg, DWI_1 ... DWI_N], all shells
bvect_all = bvect(dwi_inds, :);
bval_list_dwi = bval_list(dwi_inds);

% For tensor estimation: indices into dwi_data selecting only the highest shell
bval = max(bval_list);
highshell_inds = find(bval_list_dwi >= bval * 0.9); % indices within dwi_data (1-based)
bvect = bvect_all(highshell_inds, :);               % gradient dirs for highest shell only
fprintf('Tensor estimation will use b=%.0f with %d directions.\n', bval, length(highshell_inds));

% Default to false (will be overwritten if Section 1.2 is run successfully)
using_anat = false;

fprintf("Done with loading!\n");

% Keep: dti_all_unreg, dti_dims, bval, bvect (highest shell), bvect_all, bval_list_dwi, highshell_inds, using_anat, output_dir
clear series_path file_list num_files idx temp_info per_frame dti_pixel_spacing dti_slice_spacing temp_vol i currentFile vol_data info seq vec b0_inds dwi_inds b0_avg dwi_data unique_bvals

% Save the DTI volume, dimensions, gradient vectors, the b-value, and the anatomical flag
% bvect_all / bval_list_dwi / highshell_inds preserved for downstream flexibility
save(fullfile(output_dir, 'dti_all_unreg.mat'), 'dti_all_unreg', 'dti_dims', 'using_anat', '-v7.3');
save(fullfile(output_dir, 'bvect.mat'), 'bvect', 'bval', 'bvect_all', 'bval_list_dwi', 'highshell_inds');

