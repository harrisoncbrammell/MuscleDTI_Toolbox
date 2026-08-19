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

