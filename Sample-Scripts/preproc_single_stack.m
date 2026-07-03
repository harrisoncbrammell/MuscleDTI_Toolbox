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

%% 1.2 Load anatomical volumes (Dixon water + fat)
% Water (t1_vibe_dixon_tra_W): registration reference — muscle bright, matches b=0 DWI contrast.
% Fat   (t1_vibe_dixon_tra_F): masking only — combined with water as pseudo-PD in Section 2.2.
% Both are resized to DTI matrix dimensions here so downstream sections get consistent arrays.

[anat_file, anat_dir] = uigetfile('*.dcm', 'Select Dixon WATER image (t1_vibe_dixon_tra_W)');

if isequal(anat_file, 0)
    fprintf("No anatomical volume selected. Registration will use b=0 as reference.\n");
else
    fprintf("Loading Dixon water volume...\n");
    full_anat_path = fullfile(anat_dir, anat_file);
    anat_vol  = squeeze(double(dicomread(full_anat_path)));
    anat_info = dicominfo(full_anat_path);

    try
        anat_per_frame     = anat_info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1;
        anat_pixel_spacing = double(anat_per_frame.PixelSpacing(1));
        anat_slice_thick   = double(anat_per_frame.SliceThickness);
    catch
        anat_pixel_spacing = 0.5;
        anat_slice_thick   = 4.0;
        fprintf('WARNING: Could not read anatomical per-frame geometry. Using fallback values.\n');
    end
    fprintf('Dixon Water: in-plane = %.4f mm, slice thickness = %.1f mm, matrix = %d x %d x %d\n', ...
        anat_pixel_spacing, anat_slice_thick, anat_info.Rows, anat_info.Columns, int32(anat_info.NumberOfFrames));

    [target_rows, target_cols, target_slices, ~] = size(dti_all_unreg);
    fprintf('Resizing Dixon water to DTI matrix (%d x %d x %d)...\n', target_rows, target_cols, target_slices);
    anat_vol = imresize3(anat_vol, [target_rows, target_cols, target_slices], 'linear');

    % Load Dixon fat for pseudo-PD masking (same geometry as water, no need to re-read header)
    [fat_file, fat_dir] = uigetfile('*.dcm', 'Select Dixon FAT image (t1_vibe_dixon_tra_F)');
    if isequal(fat_file, 0)
        fprintf("No fat image selected. Section 2.2 will use Dixon water only for masking.\n");
    else
        fprintf("Loading Dixon fat volume...\n");
        fat_vol = squeeze(double(dicomread(fullfile(fat_dir, fat_file))));
        fat_vol = imresize3(fat_vol, [target_rows, target_cols, target_slices], 'linear');
        fprintf("Dixon fat loaded and resized.\n");
    end

    % After resize, anat_dims matches DTI spatial matrix
    anat_dims  = dti_dims;
    using_anat = true;
    fprintf('Anatomical volumes ready.\n');
end

clear anat_file anat_dir fat_file fat_dir target_rows target_cols target_slices full_anat_path anat_info anat_per_frame anat_pixel_spacing anat_slice_thick

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

% Use pseudo-PD (fat + water) as the masking reference, matching the
% toolbox pre_process.m approach. Pseudo-PD makes the full leg cross-section
% bright (both fat and muscle visible), giving a clean outer boundary for
% Otsu thresholding before erosion strips the skin/fat layer. Dixon water
% alone risks a jagged boundary where subcutaneous fat is dark.
% Fall back to b=0 if no anatomical volumes are loaded.
if exist('anat_vol', 'var') && exist('fat_vol', 'var')
    fprintf("Using pseudo-PD (fat + water) for masking reference.\n");
    ref_vol = anat_vol + fat_vol;
elseif exist('anat_vol', 'var')
    fprintf("Dixon fat not loaded. Using Dixon water only for masking reference.\n");
    ref_vol = anat_vol;
else
    fprintf("No anatomical volume found. Using average b=0 image for masking reference.\n");
    ref_vol = dti_all_unreg(:, :, :, 1);
end

[rows, cols, slices] = size(ref_vol);
pd_mask = zeros(rows, cols, slices);

for z = 1:slices
    loop_image = ref_vol(:,:,z);

    if max(loop_image(:)) > 0
        loop_image = loop_image / max(loop_image(:));
    end

    % Otsu threshold
    loop_mask = zeros(rows, cols);
    threshold = graythresh(loop_image);
    loop_mask(loop_image > threshold) = 1;

    % Morphological cleanup:
    % 1. Remove small noise blobs
    % 2. Fill holes selectively (invert trick avoids filling large background regions)
    % 3. Erode twice to strip skin and subcutaneous fat pixels at the boundary
    loop_mask = bwareaopen(logical(loop_mask), 500);
    inverted_mask = ~loop_mask;
    inverted_mask = bwareaopen(inverted_mask, 200);
    loop_mask = ~inverted_mask;
    loop_mask = bwmorph(loop_mask, 'erode', 2);

    pd_mask(:,:,z) = loop_mask;
end

clear z loop_image threshold loop_mask ref_vol rows cols slices inverted_mask

% Save with descriptive name following {descriptor}_trial_{N}_{method} convention
if ~exist('output_dir', 'var')
    output_dir = uigetdir(pwd, 'Select output folder to save mask');
end
save(fullfile(output_dir, 'pd_mask_trial_1_otsu.mat'), 'pd_mask', '-v7.3');
fprintf("Mask successfully generated and saved to %s\n", output_dir);

%% 2.2b Quadrant-Adaptive Otsu masking (improved edge-slice version)
% Calls form_dwi_mask (Local Functions) on each slice of the anatomical
% reference. Fixes coil-dropout holes on slices 1-4 by thresholding each
% image quadrant independently. Compare against 2.2 using Section 2.7.

fprintf('Generating mask with quadrant-adaptive Otsu (form_dwi_mask)...\n');

if exist('anat_vol', 'var') && exist('fat_vol', 'var')
    ref_vol = anat_vol + fat_vol;
elseif exist('anat_vol', 'var')
    ref_vol = anat_vol;
else
    ref_vol = dti_all_unreg(:,:,:,1);
end

n_slices = size(ref_vol, 3);
pd_mask  = false(size(ref_vol));
for z = 1:n_slices
    pd_mask(:,:,z) = form_dwi_mask(ref_vol(:,:,z));
end

clear ref_vol n_slices z

if ~exist('output_dir', 'var')
    output_dir = uigetdir(pwd, 'Select output folder to save mask');
end
save(fullfile(output_dir, 'pd_mask_trial_2_otsu_quadrant.mat'), 'pd_mask', '-v7.3');
fprintf('Quadrant-adaptive mask saved to %s\n', output_dir);

%% 2.3 Upload mask from mat file
start_dir = '';
if exist('output_dir', 'var'), start_dir = output_dir; end
[mask_file, mask_dir] = uigetfile(fullfile(start_dir, '*.mat'), 'Select the mask .mat file');

if isequal(mask_file, 0)
    fprintf("Mask load cancelled.\n");
else
    pd_mask = load(fullfile(mask_dir, mask_file));
    if isstruct(pd_mask)
        fields = fieldnames(pd_mask);
        pd_mask = pd_mask.(fields{1});
    end
    fprintf("Mask successfully loaded from %s\n", fullfile(mask_dir, mask_file));
end

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
    
    if ~exist('output_dir', 'var')
        output_dir = uigetdir(pwd, 'Select output folder to save mask');
    end
    save(fullfile(output_dir, 'pd_mask_trial_1_define_muscle.mat'), 'pd_mask', '-v7.3');
    fprintf("define_muscle complete! Mask saved to %s\n", output_dir);
end

clear total_slices prompt dlgtitle dims definput answer slices_to_segment fv_options ref_vol current_dims

%% 2.6 Mask Quality Control
% Run after any masking section (2.1-2.5) to evaluate the current pd_mask.
% Works standalone — only requires pd_mask in the workspace.

if ~exist('pd_mask', 'var')
    fprintf('No pd_mask in workspace. Run a masking section first.\n');
else

[qc_rows, qc_cols, qc_slices] = size(pd_mask);

% --- Choose reference image for overlay ---
if exist('anat_vol', 'var')
    qc_ref = anat_vol;   qc_ref_label = 'Dixon Water';
elseif exist('dti_all_unreg', 'var')
    qc_ref = dti_all_unreg(:,:,:,1);   qc_ref_label = 'Average b=0';
else
    qc_ref = [];   qc_ref_label = 'none';
end

% --- Analytical summary ---
slice_areas   = squeeze(sum(sum(pd_mask, 1), 2));  % masked voxels per slice
total_voxels  = sum(pd_mask(:));
empty_slices  = find(slice_areas == 0);
low_threshold = mean(slice_areas(slice_areas > 0)) - 2*std(slice_areas(slice_areas > 0));
outlier_slices = find(slice_areas > 0 & slice_areas < low_threshold);

% Per-slice solidity: ratio of mask area to convex hull area.
% Close to 1 = compact/convex muscle shape. Low values = ragged or fragmented mask.
solidity = zeros(qc_slices, 1);
for s = 1:qc_slices
    props = regionprops(logical(pd_mask(:,:,s)), 'Solidity');
    if ~isempty(props)
        solidity(s) = max([props.Solidity]);
    end
end

% Physical units (requires dti_dims from Section 1.1)
if exist('dti_dims', 'var')
    vox_area_mm2  = dti_dims(1)^2;
    vox_vol_mm3   = dti_dims(1)^2 * dti_dims(2);
    phys_areas    = slice_areas * vox_area_mm2;
    total_vol_mL  = total_voxels * vox_vol_mm3 / 1000;
    area_unit     = 'mm²';
else
    phys_areas   = slice_areas;
    total_vol_mL = NaN;
    area_unit    = 'voxels';
end

fprintf('=== Mask QC Summary ===\n');
fprintf('  Mask matrix:          %d x %d x %d\n', qc_rows, qc_cols, qc_slices);
fprintf('  Total masked voxels:  %d\n', total_voxels);
if ~isnan(total_vol_mL)
    fprintf('  Estimated volume:     %.1f mL\n', total_vol_mL);
end
fprintf('  Mean area per slice:  %.0f %s\n', mean(phys_areas), area_unit);
fprintf('  Mean solidity:        %.3f (1.0 = perfectly convex)\n', mean(solidity(solidity > 0)));
if ~isempty(empty_slices)
    fprintf('  WARNING - Empty slices (%d): %s\n', length(empty_slices), num2str(empty_slices'));
else
    fprintf('  Empty slices:         none\n');
end
if ~isempty(outlier_slices)
    fprintf('  WARNING - Low-area outliers >2SD below mean (%d): %s\n', length(outlier_slices), num2str(outlier_slices'));
else
    fprintf('  Low-area outliers:    none\n');
end

% --- Per-slice area + solidity figure ---
figure('Name', 'Mask QC: Per-Slice Metrics', 'Color', 'w', 'Position', [100 100 900 400]);

subplot(1,2,1);
plot(1:qc_slices, phys_areas, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.2); hold on;
yline(mean(phys_areas), 'r--', 'Mean', 'LineWidth', 1);
if ~isempty(empty_slices)
    scatter(empty_slices, zeros(size(empty_slices)), 60, 'r', 'filled');
end
if ~isempty(outlier_slices)
    scatter(outlier_slices, phys_areas(outlier_slices), 60, [1 0.5 0], 'filled');
end
xlabel('Slice'); ylabel(['Area (' area_unit ')']);
title('Masked area per slice'); grid on;

subplot(1,2,2);
plot(1:qc_slices, solidity, 'g-o', 'MarkerSize', 4, 'LineWidth', 1.2); hold on;
yline(0.9, 'r--', 'Good threshold', 'LineWidth', 1);
xlabel('Slice'); ylabel('Solidity');
title('Mask solidity per slice (1 = convex)'); ylim([0 1.05]); grid on;

% --- Visual overlay (scrollable) ---
if ~isempty(qc_ref)
    fprintf('Building overlay viewer (%s + mask)...\n', qc_ref_label);
    qc_ref_norm    = mat2gray(qc_ref);
    overlay_stack  = zeros(qc_rows, qc_cols, 3, qc_slices, 'uint8');
    for s = 1:qc_slices
        frame = labeloverlay(qc_ref_norm(:,:,s), logical(pd_mask(:,:,s)), ...
                             'Colormap', [1 0 0], 'Transparency', 0.5);
        overlay_stack(:,:,:,s) = im2uint8(frame);
    end
    implay(overlay_stack); % already [rows, cols, 3, slices] — implay RGB format
end

clear qc_rows qc_cols qc_slices qc_ref qc_ref_label qc_ref_norm overlay_stack ...
      slice_areas total_voxels empty_slices low_threshold outlier_slices ...
      solidity s props phys_areas total_vol_mL vox_area_mm2 vox_vol_mm3 area_unit frame

end % pd_mask existence check

%% 2.7 Compare two masks
% Load a second mask and compute Dice similarity — useful for comparing
% masking methods (e.g. Otsu vs define_muscle) or evaluating parameter changes.

start_dir = '';
if exist('output_dir', 'var'), start_dir = output_dir; end

[f1, d1] = uigetfile(fullfile(start_dir, '*.mat'), 'Select FIRST mask');
[f2, d2] = uigetfile(fullfile(start_dir, '*.mat'), 'Select SECOND mask');

if isequal(f1, 0) || isequal(f2, 0)
    fprintf('Comparison cancelled.\n');
else
    m1 = load(fullfile(d1, f1)); if isstruct(m1), fld = fieldnames(m1); m1 = m1.(fld{1}); end
    m2 = load(fullfile(d2, f2)); if isstruct(m2), fld = fieldnames(m2); m2 = m2.(fld{1}); end

    n_slices = size(m1, 3);
    slice_dice = zeros(n_slices, 1);
    for s = 1:n_slices
        a = logical(m1(:,:,s));  b = logical(m2(:,:,s));
        slice_dice(s) = 2*sum(a(:) & b(:)) / (sum(a(:)) + sum(b(:)) + eps);
    end
    overall_dice = 2*sum(m1(:) & m2(:)) / (sum(m1(:)) + sum(m2(:)));

    fprintf('=== Mask Comparison ===\n');
    fprintf('  Mask 1: %s\n', f1);
    fprintf('  Mask 2: %s\n', f2);
    fprintf('  Overall Dice:      %.4f\n', overall_dice);
    fprintf('  Mean per-slice Dice: %.4f\n', mean(slice_dice));
    fprintf('  Min per-slice Dice:  %.4f (slice %d)\n', min(slice_dice), find(slice_dice==min(slice_dice),1));

    figure('Name', 'Mask Comparison: Per-Slice Dice', 'Color', 'w');
    plot(1:n_slices, slice_dice, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.2); hold on;
    yline(overall_dice, 'r--', sprintf('Overall Dice = %.3f', overall_dice), 'LineWidth', 1);
    xlabel('Slice'); ylabel('Dice coefficient');
    title(sprintf('Per-slice Dice: %s vs %s', f1, f2));
    ylim([0 1.05]); grid on;

    clear m1 m2 f1 f2 d1 d2 fld n_slices s a b slice_dice overall_dice start_dir
end

%% 3.1 Demons registration method

fprintf('Starting Registration (this may take time)....\n');

% 1. Ensure the mask is a clean matrix
if isstruct(pd_mask) 
    fields = fieldnames(pd_mask);
    pd_mask = pd_mask.(fields{1});
    fprintf("Fixed pd_mask struct. Now it is a matrix.\n");
end

% setup volumes and pre-allocate
[~, ~, slcs, total_vols] = size(dti_all_unreg);
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

        % Fixed side: mask with the static, pre-computed whole-leg mask.
        % That boundary is trusted -- the anatomical/reference image isn't
        % the one being warped.
        fixed_masked = fixed_slice .* current_mask;

        % Moving side: mask with a FRESH mask built from THIS slice's own
        % contrast, not the shared pd_mask. Reusing pd_mask here would crop
        % the moving image to a boundary computed for a different image,
        % silently assuming the moving tissue edge is already where pd_mask
        % says it should be -- which is exactly what registration is
        % supposed to figure out. Quadrant-adaptive Otsu + dilate gives the
        % moving edge some margin so demons has real signal to align there.
        % See form_dwi_mask (Local Functions, bottom of script); adapted
        % from bdamon/MuscleDTI_Toolbox Sample-Scripts/pre_process.m.
        moving_mask = form_dwi_mask(moving_slice);
        moving_masked = moving_slice .* moving_mask;

        % Safely normalize (Skip registration for this slice if either mask is empty)
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

save(fullfile(output_dir, 'dti_registered_trial_1.mat'), 'dti_all_reg', 'pd_mask', 'bvect', '-v7.3');
fprintf('Registration Complete! Saved to %s\n', output_dir);

clear fields rows cols slcs total_vols fixed_vol v z fixed_slice moving_slice current_mask moving_mask fixed_masked moving_masked fixed_max moving_max fixed_norm moving_norm disp_field

%% 3.1.1 Demons registration method cpu accelorated
fprintf('Starting Registration (this may take time)....\n');

% ensure the mask is a clean matrix
if isstruct(pd_mask) 
    fields = fieldnames(pd_mask);
    pd_mask = pd_mask.(fields{1});
    fprintf("Fixed pd_mask struct. Now it is a matrix.\n");
end

% setup volumes and pre-allocate
[rows, cols, slcs, total_vols] = size(dti_all_unreg);
dti_all_reg = zeros(size(dti_all_unreg));

% pick ground truth reference image
if exist('using_anat', 'var') && using_anat && exist('anat_vol', 'var')
    fprintf('Using Anatomical volume as fixed reference for registration.\n');
    fixed_vol = anat_vol;
else
    fprintf('Using average b=0 volume as fixed reference for registration.\n');
    fixed_vol = dti_all_unreg(:, :, :, 1);
    using_anat = false; % enforce false just in case
    dti_all_reg(:, :, :, 1) = dti_all_unreg(:, :, :, 1); % b=0 is perfectly aligned with itself
end

% Pre-slice fixed_vol and pd_mask into cell arrays so parfor workers receive
% only the slices they need rather than broadcasting the entire 3D volume.
fixed_slices = cell(slcs, 1);
mask_slices  = cell(slcs, 1);
for z = 1:slcs
    fixed_slices{z} = fixed_vol(:, :, z);
    mask_slices{z}  = pd_mask(:, :, z);
end

% The Parallel Registration Loop
fprintf('Starting Parallel Slice-by-Slice Demons Registration...\n');
fprintf('Booting up parallel pool (this takes a few seconds the first time)...\n');

parfor v = 1:total_vols 
    
    % Skip registering the b=0 image to itself if we aren't using an anatomical reference
    if ~using_anat && v == 1
        continue; 
    end

    fprintf('  Core processing Volume %d of %d...\n', v, total_vols);
    
    % Create a temporary 3D volume for this specific worker core to build
    temp_vol_reg = zeros(rows, cols, slcs);
    
    for z = 1:slcs
        % Use pre-sliced cell arrays to avoid broadcasting full 3D volumes
        fixed_slice  = fixed_slices{z};
        moving_slice = dti_all_unreg(:, :, z, v);
        current_mask = mask_slices{z};

        % Fixed side: mask with the static, pre-computed whole-leg mask.
        % That boundary is trusted -- the anatomical/reference image isn't
        % the one being warped.
        fixed_masked = fixed_slice .* current_mask;

        % Moving side: mask with a FRESH mask built from THIS slice's own
        % contrast, not the shared pd_mask. See the matching comment in
        % Section 3.1 and form_dwi_mask (Local Functions, bottom of script)
        % for the rationale. form_dwi_mask is self-contained (stateless
        % besides its own input), so it's safe to call inside parfor.
        moving_mask = form_dwi_mask(moving_slice);
        moving_masked = moving_slice .* moving_mask;

        % Safely normalize (Skip registration for this slice if either mask is empty)
        fixed_max = max(fixed_masked(:));
        moving_max = max(moving_masked(:));

        if fixed_max > 0 && moving_max > 0
            fixed_norm = fixed_masked / fixed_max;
            moving_norm = moving_masked / moving_max;
        else
            temp_vol_reg(:, :, z) = moving_slice; % Copy original and skip
            continue;
        end

        % THE CONTRAST HACK: Apply only when registering b=0 to Anatomical
        if using_anat && v == 1
            moving_norm(moving_norm > 0.5) = 0.5; % Flatten out bright blood vessels
            moving_norm = moving_norm * 2;        % Scale up to match anatomical brightness
        end

        % Calculate displacement field and apply the warp to the UNMASKED original slice
        [disp_field, ~] = imregdemons(moving_norm, fixed_norm, [500 400 200], 'AccumulatedFieldSmoothing', 1.3, 'DisplayWaitbar', false);
        temp_vol_reg(:, :, z) = imwarp(moving_slice, disp_field);
    end

    % Drop the finished volume into the master array
    dti_all_reg(:, :, :, v) = temp_vol_reg;
end

save(fullfile(output_dir, 'dti_registered_trial_1.mat'), 'dti_all_reg', 'pd_mask', 'bvect', '-v7.3');
fprintf('Registration Complete! Saved to %s\n', output_dir);

clear fields rows cols slcs total_vols fixed_vol v z fixed_slice moving_slice current_mask moving_mask fixed_masked moving_masked fixed_max moving_max fixed_norm moving_norm disp_field temp_vol_reg
%% 3.2 FSL Eddy based registraton method

%TODO: base of off exhisting implementation in eddy_correct.m

%% 3.3 Registration Quality Control (Visual Inspection)

% 1. Settings: Choose which volume and slice to inspect
% (Volume 2 is usually the first diffusion direction)
v_inspect = 1; % b=0 matches Dixon water contrast best — most diagnostic; change to 2+ for a DWI direction
z_inspect = round(size(dti_all_reg, 3) / 2); % Middle slice

% 2. Prepare the images
% Grab the fixed reference used in Section 3.1
if exist('using_anat', 'var') && using_anat
    fixed_img = anat_vol(:,:,z_inspect);
    ref_label = 'Anatomical';
else
    fixed_img = dti_all_unreg(:,:,z_inspect, 1);
    ref_label = 'Average b=0';
end

moving_unreg = dti_all_unreg(:,:,z_inspect, v_inspect);
moving_reg = dti_all_reg(:,:,z_inspect, v_inspect);

% 3. Create the QC Figure
figure('Name', 'Registration Quality Control', 'NumberTitle', 'off', 'Color', 'w');

% Subplot 1: Checkerboard (Alignment Test)
subplot(2,2,1);
imshowpair(fixed_img, moving_reg, 'checkerboard');
title(['Checkerboard: ', ref_label, ' vs Reg Vol ', num2str(v_inspect)]);
xlabel('Discontinuities in lines = Bad Alignment');

% Subplot 2: Difference Map (Error Heatmap)
% We normalize both to [0 1] before subtracting to see pure alignment errors
f_norm = fixed_img / max(fixed_img(:));
m_norm = moving_reg / max(moving_reg(:));
diff_map = abs(f_norm - m_norm);

subplot(2,2,2);
imshow(diff_map, []);
colormap(gca, 'hot'); 
colorbar;
title('Difference Map (Normalized)');
xlabel('Bright outlines = Registration Error');

% Subplot 3: Overlay (Green/Magenta)
% Fixed is Green, Moving is Magenta. Grey areas are perfectly aligned.
subplot(2,2,3);
imshowpair(fixed_img, moving_reg, 'falsecolor');
title('Overlay: Green(Fixed) / Mag(Reg)');
xlabel('Color fringes indicate mismatch');

% Subplot 4: Before vs After (Side-by-Side)
subplot(2,2,4);
imshowpair(moving_unreg, moving_reg, 'montage');
title('Before (Left) vs After (Right)');
xlabel('Compare distortion/warping');

fprintf('QC Figure generated for Slice %d, Volume %d.\n', z_inspect, v_inspect);

% 4. Optional: Animation (The Flicker Test)
fprintf('Starting Flicker Test in a new window...\n');
figure('Name', 'Flicker Test: Toggle between Fixed and Registered', 'Color', 'k');
for i = 1:6
    imshow(fixed_img, []); title(['FIXED (', ref_label, ')']); pause(0.4);
    imshow(moving_reg, []); title(['REGISTERED (Vol ', num2str(v_inspect), ')']); pause(0.4);
end

%% 4.1 Denoising with anisotropic smoothing
% Applies edge-preserving anisotropic smoothing to the registered DTI volumes.
% Parameters follow Buck et al. (PLoS ONE, 2015) as used in pre_process.m.
% delta_t is scaled by the actual number of volumes in our data.
% Requires aniso4D_smoothing.m on the MATLAB path (Preprocessing-Functions/).

noise      = 5;
sigma      = noise / 100;
rho        = 2 * sigma;
delta_t    = noise * 3 / size(dti_all_reg, 4);
schemetype = 'implicit_multisplitting';
isnormg    = false;
isfasteig  = true;

% aniso4D_smoothing expects [x_spacing, y_spacing, slice_spacing] (3 elements).
% dti_dims is [in-plane pixel size, slice spacing] — expand to 3-element vector.
% In-plane pixels are isotropic so x and y spacing are the same.
dti_res = [dti_dims(1), dti_dims(1), dti_dims(2)];

fprintf('Running anisotropic smoothing (%d volumes, delta_t = %.4f)...\n', size(dti_all_reg, 4), delta_t);
dti_all_smooth = aniso4D_smoothing(dti_all_reg, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);

diff_max = max(abs(dti_all_smooth(:) - dti_all_reg(:)));
if diff_max < 1e-6
    fprintf('WARNING: Smoothing produced no detectable change (max diff = %.2e).\n', diff_max);
    fprintf('  -> Run: which(''aniso4D_smoothing'') to confirm the function is on the path.\n');
    fprintf('  -> Run: addpath(''tools/MuscleDTI_Toolbox/Preprocessing-Functions'') if needed.\n');
else
    fprintf('Smoothing complete. Max pixel change: %.2f\n', diff_max);
end
clear diff_max

save(fullfile(output_dir, 'dti_smooth_trial_1_aniso.mat'), 'dti_all_smooth', '-v7.3');
fprintf('Saved to %s\n', output_dir);

clear noise sigma rho delta_t schemetype isnormg isfasteig dti_res

%% 4.2 Threshold based principle component analysis (PCA) denoising

%TODO: Implement based of script provided

%% 4.3 Denoising QC — single method
% Evaluates one denoised volume against raw registered data in one figure.
% Uses dti_all_smooth if in workspace, otherwise prompts for a .mat file.
% Requires dti_all_reg in workspace.

if ~exist('dti_all_reg', 'var')
    fprintf('dti_all_reg not in workspace. Load registered data first.\n');
else

start_dir = '';
if exist('output_dir', 'var'), start_dir = output_dir; end

if ~exist('dti_all_smooth', 'var')
    [f1, d1] = uigetfile(fullfile(start_dir, '*.mat'), 'Select denoised volume');
    if isequal(f1, 0), error('No file selected.'); end
    tmp = load(fullfile(d1, f1)); flds = fieldnames(tmp);
    dti_all_smooth = tmp.(flds{1}); smooth_label = strrep(f1, '_', '\_');
    clear tmp flds f1 d1
else
    smooth_label = 'dti\_all\_smooth (workspace)';
end

z_qc = round(size(dti_all_reg, 3) / 2);
v_qc = size(dti_all_reg, 4);

% tSNR
tsnr_raw    = mean(dti_all_reg,    4) ./ (std(dti_all_reg,    0, 4) + eps);
tsnr_smooth = mean(dti_all_smooth, 4) ./ (std(dti_all_smooth, 0, 4) + eps);
tsnr_diff   = tsnr_smooth(:,:,z_qc) - tsnr_raw(:,:,z_qc);
clim_tsnr   = [0 max(prctile(tsnr_raw(:), 99), 0.01)];
clim_diff   = max(abs(tsnr_diff(:)));
if clim_diff == 0, clim_diff = 0.01; end

if exist('pd_mask', 'var')
    m = logical(pd_mask);
    r_raw = mean(tsnr_raw(m)); r_sm = mean(tsnr_smooth(m));
    fprintf('=== Denoising QC: tSNR inside mask ===\n');
    fprintf('  Raw:      %.2f\n', r_raw);
    fprintf('  Denoised: %.2f  (%+.1f%%)\n', r_sm, 100*(r_sm - r_raw)/r_raw);
    clear m r_raw r_sm
end

% Residual
residual  = dti_all_reg(:,:,z_qc,v_qc) - dti_all_smooth(:,:,z_qc,v_qc);
res_lim   = max(abs(residual(:)));
if res_lim == 0, res_lim = 1; end   % guard against all-zero residual (no smoothing applied)

% Central mask voxel for signal profile
if exist('pd_mask', 'var')
    [ri, ci] = find(logical(pd_mask(:,:,z_qc)));
    v_row = round(median(ri)); v_col = round(median(ci));
else
    v_row = round(size(dti_all_reg,1)/2); v_col = round(size(dti_all_reg,2)/2);
end
sig_raw    = squeeze(dti_all_reg(v_row,    v_col, z_qc, :));
sig_smooth = squeeze(dti_all_smooth(v_row, v_col, z_qc, :));

% Single figure: 3 rows × 3 cols
% Colorbars omitted on image panels — scale shown in title to avoid layout squishing.
figure('Name', 'Denoising QC', 'Color', 'w', 'Position', [50 50 1100 750]);

no_change = (res_lim == 1 && max(abs(residual(:))) == 0);  % true when smoothing did nothing

% Row 1: tSNR maps
ax1 = subplot(3,3,1);
imagesc(ax1, tsnr_raw(:,:,z_qc), clim_tsnr); colormap(ax1, 'hot');
set(ax1, 'Visible', 'off', 'DataAspectRatio', [1 1 1]);
title(ax1, sprintf('tSNR raw (slice %d)  [0–%.1f]', z_qc, clim_tsnr(2)), 'Color', 'k', 'Visible', 'on');

ax2 = subplot(3,3,2);
imagesc(ax2, tsnr_smooth(:,:,z_qc), clim_tsnr); colormap(ax2, 'hot');
set(ax2, 'Visible', 'off', 'DataAspectRatio', [1 1 1]);
title(ax2, sprintf('tSNR denoised — %s', smooth_label), 'Color', 'k', 'Visible', 'on');

ax3 = subplot(3,3,3);
show_or_warn(ax3, tsnr_diff, [-clim_diff clim_diff], bwr_colormap(), ...
    sprintf('tSNR diff  [±%.2f]', clim_diff), no_change);

% Row 2: image slices + residual
ax4 = subplot(3,3,4);
imagesc(ax4, dti_all_reg(:,:,z_qc,v_qc)); colormap(ax4, 'gray');
set(ax4, 'Visible', 'off', 'DataAspectRatio', [1 1 1]);
title(ax4, sprintf('Raw (slice %d, vol %d)', z_qc, v_qc), 'Color', 'k', 'Visible', 'on');

ax5 = subplot(3,3,5);
imagesc(ax5, dti_all_smooth(:,:,z_qc,v_qc)); colormap(ax5, 'gray');
set(ax5, 'Visible', 'off', 'DataAspectRatio', [1 1 1]);
title(ax5, 'Denoised', 'Color', 'k', 'Visible', 'on');

ax6 = subplot(3,3,6);
show_or_warn(ax6, residual, [-res_lim res_lim], bwr_colormap(), ...
    'Residual — no tissue edges = good', no_change);

% Row 3: signal profile spanning full width
ax7 = subplot(3,3,7:9);
plot(ax7, sig_raw,    'b-', 'LineWidth', 1,   'DisplayName', 'Raw'); hold(ax7, 'on');
plot(ax7, sig_smooth, 'r-', 'LineWidth', 1.5, 'DisplayName', smooth_label);
if no_change
    text(ax7, mean(xlim(ax7)), mean(ylim(ax7)), 'WARNING: Raw = Denoised — smoothing had no effect', ...
        'HorizontalAlignment', 'center', 'Color', [0.6 0 0], 'FontSize', 10);
end
xlabel(ax7, 'Volume index', 'Color', 'k');
ylabel(ax7, 'Signal intensity', 'Color', 'k');
title(ax7, sprintf('Signal profile — voxel [%d,%d] slice %d  |  smoother = good, same shape = no over-smoothing', v_row, v_col, z_qc), 'Color', 'k');
set(ax7, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
legend(ax7); grid(ax7, 'on');

clear tsnr_raw tsnr_smooth tsnr_diff clim_tsnr clim_diff residual res_lim ...
      sig_raw sig_smooth v_row v_col z_qc v_qc ri ci smooth_label start_dir ...
      ax1 ax2 ax3 ax4 ax5 ax6 ax7 no_change

end

%% 4.4 Denoising QC — compare two methods
% Loads two denoised .mat files and compares them against each other and raw.
% Requires dti_all_reg in workspace.

if ~exist('dti_all_reg', 'var')
    fprintf('dti_all_reg not in workspace. Load registered data first.\n');
else

start_dir = '';
if exist('output_dir', 'var'), start_dir = output_dir; end

[fa, da] = uigetfile(fullfile(start_dir, '*.mat'), 'Select method A denoised volume');
[fb, db] = uigetfile(fullfile(start_dir, '*.mat'), 'Select method B denoised volume');
if isequal(fa, 0) || isequal(fb, 0), error('Two files required for comparison.'); end

tmp = load(fullfile(da, fa)); flds = fieldnames(tmp); smooth_a = tmp.(flds{1}); clear tmp flds;
tmp = load(fullfile(db, fb)); flds = fieldnames(tmp); smooth_b = tmp.(flds{1}); clear tmp flds;

z_qc  = round(size(dti_all_reg, 3) / 2);
v_qc  = size(dti_all_reg, 4);

tsnr_raw = mean(dti_all_reg, 4) ./ (std(dti_all_reg, 0, 4) + eps);
tsnr_a   = mean(smooth_a,    4) ./ (std(smooth_a,    0, 4) + eps);
tsnr_b   = mean(smooth_b,    4) ./ (std(smooth_b,    0, 4) + eps);

if exist('pd_mask', 'var')
    m = logical(pd_mask);
    raw_mean = mean(tsnr_raw(m));
    fprintf('=== Denoising Comparison: tSNR inside mask ===\n');
    fprintf('  Raw:      %.2f\n', raw_mean);
    fprintf('  Method A (%s): %.2f  (%.1f%% improvement)\n', fa, mean(tsnr_a(m)), 100*(mean(tsnr_a(m))-raw_mean)/raw_mean);
    fprintf('  Method B (%s): %.2f  (%.1f%% improvement)\n', fb, mean(tsnr_b(m)), 100*(mean(tsnr_b(m))-raw_mean)/raw_mean);
    clear m raw_mean
end

clim = [0 prctile(tsnr_raw(:), 99)];
figure('Name', 'Denoising Comparison: tSNR', 'Color', 'w', 'Position', [50 50 1400 400]);
subplot(1,4,1); imagesc(tsnr_raw(:,:,z_qc), clim); colormap(gca,'hot'); colorbar;
title(sprintf('Raw (slice %d)', z_qc)); axis image off;
subplot(1,4,2); imagesc(tsnr_a(:,:,z_qc), clim); colormap(gca,'hot'); colorbar;
title(sprintf('Method A — %s', fa)); axis image off;
subplot(1,4,3); imagesc(tsnr_b(:,:,z_qc), clim); colormap(gca,'hot'); colorbar;
title(sprintf('Method B — %s', fb)); axis image off;
subplot(1,4,4); imagesc(tsnr_b(:,:,z_qc) - tsnr_a(:,:,z_qc)); colormap(gca, bwr_colormap()); colorbar;
title('Difference (B − A)'); axis image off;

% Signal profile comparison
if exist('pd_mask', 'var')
    [ri, ci] = find(logical(pd_mask(:,:,z_qc)));
    v_row = round(median(ri)); v_col = round(median(ci));
else
    v_row = round(size(dti_all_reg,1)/2); v_col = round(size(dti_all_reg,2)/2);
end

figure('Name', 'Denoising Comparison: Signal Profile', 'Color', 'w', 'Position', [50 100 800 350]);
plot(squeeze(dti_all_reg(v_row,v_col,z_qc,:)), 'b-',  'LineWidth', 1,   'DisplayName', 'Raw'); hold on;
plot(squeeze(smooth_a(v_row,v_col,z_qc,:)),    'r-',  'LineWidth', 1.5, 'DisplayName', fa);
plot(squeeze(smooth_b(v_row,v_col,z_qc,:)),    'g-',  'LineWidth', 1.5, 'DisplayName', fb);
xlabel('Volume index'); ylabel('Signal intensity');
title(sprintf('Signal profile — voxel [%d, %d, slice %d]', v_row, v_col, z_qc));
legend; grid on;

% Residual comparison
res_a = dti_all_reg(:,:,z_qc,v_qc) - smooth_a(:,:,z_qc,v_qc);
res_b = dti_all_reg(:,:,z_qc,v_qc) - smooth_b(:,:,z_qc,v_qc);
clim_r = [-1 1] * max(max(abs(res_a(:))), max(abs(res_b(:))));
figure('Name', 'Denoising Comparison: Residuals', 'Color', 'w', 'Position', [50 500 800 350]);
subplot(1,2,1); imagesc(res_a, clim_r); colormap(gca, bwr_colormap()); colorbar;
title(sprintf('Residual — %s', fa)); axis image off;
subplot(1,2,2); imagesc(res_b, clim_r); colormap(gca, bwr_colormap()); colorbar;
title(sprintf('Residual — %s', fb)); axis image off;

clear tsnr_raw tsnr_a tsnr_b smooth_a smooth_b clim clim_r ...
      res_a res_b v_row v_col z_qc v_qc ri ci fa fb da db start_dir

end

%% 5.1 Tensor estimation — weighted least-squares fit (signal2tensor2)
% Fits the Stejskal-Tanner equation per masked voxel using the averaged b=0
% and highest-shell DWI volumes. Follows the pre_process.m toolbox example
% but restructures the loop order (slice outermost) so the inner row/col
% loops run in parallel across cores via parfor.
%
% Inputs required in workspace:
%   dti_all_smooth  (preferred) or dti_all_reg (if denoising was skipped)
%   pd_mask, bvect, bval, highshell_inds, output_dir
%
% Output:
%   tensor_m — [rows × cols × slices × 3 × 3] diffusion tensor field (units: 1/b = mm²/s)

% Use denoised data if available, otherwise fall back to registered data
if exist('dti_all_smooth', 'var')
    dti_input = dti_all_smooth;
    input_label = 'dti_all_smooth';
elseif exist('dti_all_reg', 'var')
    fprintf('WARNING: dti_all_smooth not found. Using undenoised dti_all_reg.\n');
    dti_input = dti_all_reg;
    input_label = 'dti_all_reg';
else
    [f, d] = uigetfile(fullfile(output_dir, '*.mat'), 'Select registered or denoised DTI volume');
    if isequal(f, 0), error('No file selected.'); end
    tmp = load(fullfile(d, f)); flds = fieldnames(tmp);
    dti_input = [];
    for fi = 1:numel(flds)
        if ndims(tmp.(flds{fi})) == 4
            dti_input = tmp.(flds{fi}); break;
        end
    end
    if isempty(dti_input), error('Could not find a 4D array in the selected file.'); end
    input_label = f;
    clear tmp flds f d fi
end

[rows, cols, slices, ~] = size(dti_input);

% Build index vector into dti_input for signal2tensor2.
% Vol 1 is the averaged b=0; the DWI block starts at vol 2.
% highshell_inds indexes within the DWI block (1-based), so +1 shifts
% them to the correct positions in dti_input.
vol_inds = [1; highshell_inds(:) + 1];
n_dirs   = length(highshell_inds);

fprintf('Tensor estimation from: %s\n', input_label);
fprintf('  %d slices, %d DWI directions (b=%.0f s/mm²)\n', slices, n_dirs, bval);
fprintf('  Fitting only within pd_mask (%d masked voxels total)\n', sum(pd_mask(:)));

% Pre-slice into cell arrays — avoids broadcasting the full 4D volume to
% every parfor worker; each worker only receives its own slice.
data_slices = cell(slices, 1);
mask_slices = cell(slices, 1);
for z = 1:slices
    data_slices{z} = dti_input(:, :, z, vol_inds);  % [rows × cols × (1+n_dirs)]
    mask_slices{z} = logical(pd_mask(:, :, z));
end
clear dti_input

% Tensor estimation — parallel over slices
% Each parfor worker fits all voxels within its slice independently.
fprintf('Running tensor estimation (parfor over slices)...\n');
tensor_slices = cell(slices, 1);

parfor z = 1:slices
    slice_data = data_slices{z};    % [rows × cols × (1+n_dirs)]
    slice_mask = mask_slices{z};    % [rows × cols] logical
    t_slice    = zeros(rows, cols, 3, 3);

    for r = 1:rows
        for c = 1:cols
            if slice_mask(r, c)
                signal_v = squeeze(slice_data(r, c, :));
                t_slice(r, c, :, :) = signal2tensor2(signal_v, bvect, bval);
            end
        end
    end

    tensor_slices{z} = t_slice;
    fprintf('  Slice %d/%d done.\n', z, slices);
end

% Reassemble slices into the full 5D tensor volume
tensor_m = zeros(rows, cols, slices, 3, 3);
for z = 1:slices
    tensor_m(:, :, z, :, :) = tensor_slices{z};
end

fprintf('Tensor estimation complete.\n');

save(fullfile(output_dir, 'tensor_trial_1.mat'), 'tensor_m', 'bvect', 'bval', 'pd_mask', '-v7.3');
fprintf('Saved tensor_trial_1.mat to %s\n', output_dir);

clear data_slices mask_slices tensor_slices t_slice slice_data slice_mask ...
      signal_v vol_inds n_dirs rows cols slices z r c input_label

%% 5.2 Tensor-derived maps and QC figure
% Computes eigenvalues/eigenvectors per masked voxel and derives FA, MD,
% TRACE, color FA, and a non-physical tensor flag. Saves all maps to a
% .mat file and displays a summary figure.
%
% Inputs:  tensor_m, pd_mask, output_dir
% Outputs: FA_map, ADC_map, trace_map, lambda1/2/3, V1_map, cFA_map, neg_eig

if ~exist('tensor_m', 'var')
    [f, d] = uigetfile(fullfile(output_dir, '*.mat'), 'Select tensor_trial_1.mat');
    if isequal(f, 0), error('tensor_m is required.'); end
    tmp = load(fullfile(d, f));
    tensor_m = tmp.tensor_m;
    clear tmp f d
end

[rows, cols, slices, ~, ~] = size(tensor_m);

% Pre-slice for parfor
t_slices = cell(slices, 1);
m_slices = cell(slices, 1);
for z = 1:slices
    t_slices{z} = squeeze(tensor_m(:, :, z, :, :));  % [rows × cols × 3 × 3]
    m_slices{z} = logical(pd_mask(:, :, z));
end

% Per-slice output cell arrays (avoids parfor classification issues)
fa_sl  = cell(slices, 1);  adc_sl  = cell(slices, 1);
tr_sl  = cell(slices, 1);  l1_sl  = cell(slices, 1);
l2_sl  = cell(slices, 1);  l3_sl  = cell(slices, 1);
v1_sl  = cell(slices, 1);  ng_sl  = cell(slices, 1);

fprintf('Computing eigendecomposition (%d slices, parfor)...\n', slices);

parfor z = 1:slices
    T  = t_slices{z};
    mk = m_slices{z};

    fa_s = zeros(rows, cols);  md_s = zeros(rows, cols);
    tr_s = zeros(rows, cols);  l1_s = zeros(rows, cols);
    l2_s = zeros(rows, cols);  l3_s = zeros(rows, cols);
    v1_s = zeros(rows, cols, 3);
    ng_s = false(rows, cols);

    for r = 1:rows
        for c = 1:cols
            if mk(r, c)
                D = squeeze(T(r, c, :, :));
                [V, L]       = eig(D);
                lam          = diag(L);
                [lam, idx]   = sort(lam, 'descend');
                V            = V(:, idx);

                md  = mean(lam);
                fa  = sqrt(3/2) * sqrt(sum((lam - md).^2)) / (sqrt(sum(lam.^2)) + eps);
                fa  = min(max(fa, 0), 1);

                fa_s(r,c)   = fa;
                md_s(r,c)   = md;
                tr_s(r,c)   = sum(lam);
                l1_s(r,c)   = lam(1);
                l2_s(r,c)   = lam(2);
                l3_s(r,c)   = lam(3);
                v1_s(r,c,:) = V(:,1);
                ng_s(r,c)   = any(lam < 0);
            end
        end
    end

    fa_sl{z} = fa_s;  adc_sl{z} = md_s;  tr_sl{z} = tr_s;
    l1_sl{z} = l1_s;  l2_sl{z} = l2_s;  l3_sl{z} = l3_s;
    v1_sl{z} = v1_s;  ng_sl{z} = ng_s;
    fprintf('  Slice %d/%d done.\n', z, slices);
end

% Reassemble into full volumes
FA_map    = zeros(rows, cols, slices);
ADC_map    = zeros(rows, cols, slices);
trace_map = zeros(rows, cols, slices);
lambda1   = zeros(rows, cols, slices);
lambda2   = zeros(rows, cols, slices);
lambda3   = zeros(rows, cols, slices);
V1_map    = zeros(rows, cols, slices, 3);
neg_eig   = false(rows, cols, slices);

for z = 1:slices
    FA_map(:,:,z)    = fa_sl{z};   ADC_map(:,:,z)    = adc_sl{z};
    trace_map(:,:,z) = tr_sl{z};   lambda1(:,:,z)   = l1_sl{z};
    lambda2(:,:,z)   = l2_sl{z};   lambda3(:,:,z)   = l3_sl{z};
    V1_map(:,:,z,:)  = v1_sl{z};   neg_eig(:,:,z)   = ng_sl{z};
end

% Color FA: FA × |first eigenvector|, channels = [R=L/R, G=A/P, B=S/I]
cFA_map = zeros(rows, cols, slices, 3);
for ch = 1:3
    cFA_map(:,:,:,ch) = FA_map .* abs(V1_map(:,:,:,ch));
end

clear t_slices m_slices fa_sl adc_sl tr_sl l1_sl l2_sl l3_sl v1_sl ng_sl ...
      z r c D V L lam idx md fa T mk ...
      fa_s md_s tr_s l1_s l2_s l3_s v1_s ng_s

% --- Command-window summary ---
mask_l = logical(pd_mask);
n_mask  = sum(mask_l(:));
n_neg   = sum(neg_eig(:));
fprintf('=== Tensor Map QC Summary ===\n');
fprintf('  Masked voxels:         %d\n',     n_mask);
fprintf('  Non-physical tensors:  %d (%.1f%%)\n', n_neg, 100*n_neg/n_mask);
fprintf('  Mean FA (±SD):         %.3f ± %.3f\n',  mean(FA_map(mask_l)), std(FA_map(mask_l)));
fprintf('  Mean ADC:              %.4f mm²/s\n',   mean(ADC_map(mask_l)));
fprintf('  Mean TRACE:            %.4f mm²/s\n',   mean(trace_map(mask_l)));
fprintf('  Expected for muscle:   FA 0.2–0.4, ADC 1.5–2.0×10⁻³, TRACE 4.5–6.0×10⁻³ mm²/s\n');

% Per-slice mean FA ± SD
slice_fa_mean = zeros(slices, 1);
slice_fa_std  = zeros(slices, 1);
for z = 1:slices
    mk_z = mask_l(:,:,z);
    if any(mk_z(:))
        fa_z = FA_map(:,:,z);
        slice_fa_mean(z) = mean(fa_z(mk_z));
        slice_fa_std(z)  = std(fa_z(mk_z));
    end
end

% --- QC figure: 2 rows × 4 cols ---
z_mid = round(slices / 2);

figure('Name', 'Tensor QC: Derived Maps', 'Color', 'w', 'Position', [50 50 1400 650]);

ax1 = subplot(2,4,1);
imagesc(ax1, FA_map(:,:,z_mid), [0 1]); colormap(ax1, 'hot'); colorbar(ax1);
set(ax1, 'DataAspectRatio', [1 1 1]); axis(ax1, 'off');
title(ax1, sprintf('FA — slice %d', z_mid), 'Color', 'k');

ax2 = subplot(2,4,2);
imagesc(ax2, ADC_map(:,:,z_mid), [0 3e-3]); colormap(ax2, 'parula'); colorbar(ax2);
set(ax2, 'DataAspectRatio', [1 1 1]); axis(ax2, 'off');
title(ax2, 'ADC = TRACE/3 (mm²/s)  [0–3×10⁻³]', 'Color', 'k');

ax3 = subplot(2,4,3);
cfa_slice = squeeze(cFA_map(:,:,z_mid,:));  % [rows × cols × 3]
image(ax3, cfa_slice);
set(ax3, 'DataAspectRatio', [1 1 1]); axis(ax3, 'off');
title(ax3, 'Color FA  (R=L/R  G=A/P  B=S/I)', 'Color', 'k');

ax4 = subplot(2,4,4);
neg_overlay = labeloverlay(mat2gray(FA_map(:,:,z_mid)), neg_eig(:,:,z_mid), ...
    'Colormap', [1 0 0], 'Transparency', 0.3);
image(ax4, neg_overlay);
set(ax4, 'DataAspectRatio', [1 1 1]); axis(ax4, 'off');
title(ax4, sprintf('Non-physical voxels (red) — %.1f%%', 100*n_neg/n_mask), 'Color', 'k');

ax5 = subplot(2,4,5);
imagesc(ax5, lambda1(:,:,z_mid)); colormap(ax5, 'parula'); colorbar(ax5);
set(ax5, 'DataAspectRatio', [1 1 1]); axis(ax5, 'off');
title(ax5, 'λ₁  (primary diffusivity)', 'Color', 'k');

ax6 = subplot(2,4,6);
imagesc(ax6, lambda2(:,:,z_mid)); colormap(ax6, 'parula'); colorbar(ax6);
set(ax6, 'DataAspectRatio', [1 1 1]); axis(ax6, 'off');
title(ax6, 'λ₂', 'Color', 'k');

ax7 = subplot(2,4,7);
imagesc(ax7, trace_map(:,:,z_mid)); colormap(ax7, 'parula'); colorbar(ax7);
set(ax7, 'DataAspectRatio', [1 1 1]); axis(ax7, 'off');
title(ax7, 'TRACE = λ₁+λ₂+λ₃', 'Color', 'k');

ax8 = subplot(2,4,8);
errorbar(ax8, 1:slices, slice_fa_mean, slice_fa_std, 'b-o', 'MarkerSize', 3, 'LineWidth', 1);
hold(ax8, 'on');
yline(ax8, 0.2, 'g--', 'Min expected', 'LineWidth', 1);
yline(ax8, 0.4, 'r--', 'Max expected', 'LineWidth', 1);
xlabel(ax8, 'Slice', 'Color', 'k'); ylabel(ax8, 'Mean FA ± SD', 'Color', 'k');
title(ax8, 'Per-slice FA profile', 'Color', 'k');
set(ax8, 'XColor', 'k', 'YColor', 'k'); grid(ax8, 'on'); ylim(ax8, [0 1]);

% --- Save all maps ---
% Save all maps — send tensor_maps_trial_1.mat to Dr. Bashir
save(fullfile(output_dir, 'tensor_maps_trial_1.mat'), ...
    'FA_map', 'ADC_map', 'trace_map', ...
    'lambda1', 'lambda2', 'lambda3', ...
    'V1_map', 'cFA_map', 'neg_eig', '-v7.3');
fprintf('Maps saved to tensor_maps_trial_1.mat in %s\n', output_dir);

clear n_mask n_neg mask_l cfa_slice neg_overlay z_mid z mk_z fa_z ...
      slice_fa_mean slice_fa_std ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8

%% 5.3 Signal residual map (tensor fit quality)
% Reconstructs the predicted DWI signal from each voxel's tensor using the
% Stejskal-Tanner forward model:
%
%   S_pred(d) = S0 * exp(-b * g_d' * D * g_d)
%
% The RMS difference between predicted and measured signals per voxel
% quantifies how well the tensor fits the data. Uniform low residuals =
% good fit. Bright patches in the residual map = bad fit (motion, signal
% dropout, or gradient non-linearity on those directions/slices).
%
% Inputs:  tensor_m, dti_all_smooth (or dti_all_reg), pd_mask, bvect, bval,
%          highshell_inds, output_dir
% Outputs: resid_map (appended to tensor_maps_trial_1.mat)

if ~exist('tensor_m', 'var')
    [f, d] = uigetfile(fullfile(output_dir, '*.mat'), 'Select tensor_trial_1.mat');
    if isequal(f, 0), error('tensor_m is required.'); end
    tmp = load(fullfile(d, f)); tensor_m = tmp.tensor_m; clear tmp f d
end

if exist('dti_all_smooth', 'var')
    dti_input = dti_all_smooth;
elseif exist('dti_all_reg', 'var')
    fprintf('WARNING: dti_all_smooth not found — using undenoised dti_all_reg.\n');
    dti_input = dti_all_reg;
else
    [f, d] = uigetfile(fullfile(output_dir, '*.mat'), 'Select registered/denoised DTI volume');
    if isequal(f, 0), error('DTI volume is required.'); end
    tmp = load(fullfile(d, f)); flds = fieldnames(tmp);
    dti_input = [];
    for fi = 1:numel(flds)
        if ndims(tmp.(flds{fi})) == 4, dti_input = tmp.(flds{fi}); break; end
    end
    if isempty(dti_input), error('No 4D array found in selected file.'); end
    clear tmp flds f d fi
end

[rows, cols, slices, ~] = size(dti_input);
vol_inds = [1; highshell_inds(:) + 1];
n_dirs   = length(highshell_inds);

% Pre-slice for parfor
data_sl   = cell(slices, 1);
tens_sl   = cell(slices, 1);
mask_sl   = cell(slices, 1);
for z = 1:slices
    data_sl{z} = dti_input(:, :, z, vol_inds);
    tens_sl{z} = squeeze(tensor_m(:, :, z, :, :));
    mask_sl{z} = logical(pd_mask(:, :, z));
end
clear dti_input

rms_sl = cell(slices, 1);
fprintf('Computing signal residuals (parfor over slices)...\n');

parfor z = 1:slices
    sig_s  = data_sl{z};   % [rows × cols × (1+n_dirs)]
    D_s    = tens_sl{z};   % [rows × cols × 3 × 3]
    mk     = mask_sl{z};
    rms_s  = zeros(rows, cols);

    for r = 1:rows
        for c = 1:cols
            if mk(r, c)
                D        = squeeze(D_s(r, c, :, :));
                signal_v = squeeze(sig_s(r, c, :));
                S0       = signal_v(1);

                S_pred = zeros(n_dirs, 1);
                for dir = 1:n_dirs
                    g          = bvect(dir, :)';
                    S_pred(dir) = S0 * exp(-bval * (g' * D * g));
                end

                S_meas     = signal_v(2:end);
                rms_s(r,c) = sqrt(mean((S_meas - S_pred).^2));
            end
        end
    end
    rms_sl{z} = rms_s;
    fprintf('  Slice %d/%d done.\n', z, slices);
end

resid_map = zeros(rows, cols, slices);
for z = 1:slices
    resid_map(:,:,z) = rms_sl{z};
end

% --- Command-window summary ---
mask_l     = logical(pd_mask);
resid_vals = resid_map(mask_l);
fprintf('=== Signal Residual QC ===\n');
fprintf('  Mean RMS residual (in mask): %.2f a.u.\n', mean(resid_vals));
fprintf('  Median RMS residual:         %.2f a.u.\n', median(resid_vals));
fprintf('  95th percentile:             %.2f a.u.\n', prctile(resid_vals, 95));
fprintf('  Max RMS residual:            %.2f a.u.\n', max(resid_vals));

% Per-slice mean residual
slice_resid = zeros(slices, 1);
for z = 1:slices
    rv = resid_map(:,:,z); mk_z = mask_l(:,:,z);
    if any(mk_z(:)), slice_resid(z) = mean(rv(mk_z)); end
end

% --- Figure: 3 panels ---
z_mid  = round(slices / 2);
clim_r = [0, prctile(resid_vals, 99)];

figure('Name', 'Tensor QC: Signal Residuals', 'Color', 'w', 'Position', [100 100 1200 380]);

ax1 = subplot(1,3,1);
imagesc(ax1, resid_map(:,:,z_mid), clim_r); colormap(ax1, 'hot'); colorbar(ax1);
set(ax1, 'DataAspectRatio', [1 1 1]); axis(ax1, 'off');
title(ax1, sprintf('RMS Residual — slice %d', z_mid), 'Color', 'k');

ax2 = subplot(1,3,2);
imagesc(ax2, FA_map(:,:,z_mid), [0 1]); colormap(ax2, 'hot'); colorbar(ax2);
set(ax2, 'DataAspectRatio', [1 1 1]); axis(ax2, 'off');
title(ax2, 'FA (same slice — for comparison)', 'Color', 'k');

ax3 = subplot(1,3,3);
plot(ax3, 1:slices, slice_resid, 'r-o', 'MarkerSize', 3, 'LineWidth', 1.2);
xlabel(ax3, 'Slice', 'Color', 'k'); ylabel(ax3, 'Mean RMS residual', 'Color', 'k');
title(ax3, 'Per-slice residual — flat curve = uniform fit quality', 'Color', 'k');
set(ax3, 'XColor', 'k', 'YColor', 'k'); grid(ax3, 'on');

% Append residual map to the existing maps file
save(fullfile(output_dir, 'tensor_maps_trial_1.mat'), 'resid_map', '-append');
fprintf('resid_map appended to tensor_maps_trial_1.mat in %s\n', output_dir);

clear data_sl tens_sl mask_sl rms_sl rms_s sig_s D_s mk D signal_v S0 S_pred S_meas ...
      g dir z r c mask_l resid_vals slice_resid rv mk_z vol_inds n_dirs ...
      rows cols slices clim_r z_mid ax1 ax2 ax3

%% Local Functions
% MATLAB requires local functions in a script file to live at the very end,
% after all executable section code -- hence this section is last regardless
% of where it's used from.

function show_or_warn(ax, data, clim, cmap, ttl, warn)
% SHOW_OR_WARN  Display an imagesc panel, or a warning text if data is unchanged.
% Used in Section 4.3 to handle the case where denoising produced no effect.
if warn
    set(ax, 'Color', [0.95 0.95 0.95], 'XColor', 'k', 'YColor', 'k');
    axis(ax, 'off');
    text(ax, 0.5, 0.5, sprintf('No change detected\n(smoothing had no effect)'), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Color', [0.6 0 0], 'FontSize', 10);
else
    imagesc(ax, data, clim); colormap(ax, cmap);
    set(ax, 'Visible', 'off', 'DataAspectRatio', [1 1 1]);
end
title(ax, ttl, 'Color', 'k', 'Visible', 'on', 'FontSize', 9);
end

function cmap = bwr_colormap()
% Blue-white-red diverging colormap (replacement for 'rdbu' which is not
% a built-in MATLAB colormap). Blue = negative, white = zero, red = positive.
n = 64;
cmap = [linspace(0,1,n)' linspace(0,1,n)' ones(n,1); ...
        ones(n,1) linspace(1,0,n)' linspace(1,0,n)'];
end

function mask = form_dwi_mask(img)
% FORM_DWI_MASK  Adaptive per-slice mask for a single moving (diffusion-
% weighted) image used in Sections 3.1 / 3.1.1 registration.
%
% Unlike pd_mask (Section 2.2), which is computed once from the anatomical/
% b0 reference and reused unchanged as the FIXED-side mask, this mask is
% meant to be recomputed fresh for every moving slice, from that slice's own
% contrast. Otsu's threshold is run independently per image quadrant (rather
% than once globally) to tolerate signal inhomogeneity across the FOV -- e.g.
% one leg sitting closer to the coil than the other, or general shading --
% then the quadrant masks are recombined, hole-filled, opened to remove
% speckle, and dilated once to give the tissue edge a small margin. That
% margin matters here specifically: it keeps the true (possibly
% EPI-distorted) moving tissue boundary inside the masked region, so demons
% still has real signal there to align, instead of a hard cutoff that hides
% the mismatch before registration even runs.
%
% Adapted from the quadrant-Otsu masking block in
% bdamon/MuscleDTI_Toolbox Sample-Scripts/pre_process.m (registration
% section), generalized here to non-square matrices and wrapped as a
% reusable function since preproc_single_stack.m calls it once per
% slice/volume rather than inline.

[rows, cols] = size(img);
half_r = floor(rows / 2);
half_c = floor(cols / 2);

% quadrant index ranges: [top-left; top-right; bottom-left; bottom-right]
idx_r = [1 half_r; 1 half_r; (half_r + 1) rows; (half_r + 1) rows];
idx_c = [1 half_c; (half_c + 1) cols; 1 half_c; (half_c + 1) cols];

quad_mask = zeros(rows, cols, 4);
for k = 1:4
    r_range = idx_r(k, 1):idx_r(k, 2);
    c_range = idx_c(k, 1):idx_c(k, 2);

    quarter_img = zeros(rows, cols);
    quarter_img(r_range, c_range) = img(r_range, c_range);

    % Guard against an all-zero quadrant (e.g. corner outside the leg)
    quarter_max = max(quarter_img(:));
    if quarter_max > 0
        quarter_img = quarter_img / quarter_max;
    end

    threshold = graythresh(quarter_img); % Otsu's method, this quadrant only
    quarter_mask = zeros(rows, cols);
    quarter_mask(quarter_img > threshold) = 1;
    quad_mask(:, :, k) = quarter_mask;
end

mask = sum(quad_mask, 3) > 0; % recombine the four quadrants
mask = imfill(mask, 'holes');
mask = bwmorph(mask, 'open');   % remove speckle noise
mask = bwmorph(mask, 'dilate'); % one pass -- give the tissue edge margin

end