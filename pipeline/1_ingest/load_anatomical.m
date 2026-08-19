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

