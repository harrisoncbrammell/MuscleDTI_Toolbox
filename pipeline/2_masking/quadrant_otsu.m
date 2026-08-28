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

% Otsu-based masking separates tissue from air, not thigh from everything
% else -- on real study data this pulls in torso and any QC phantom
% (e.g. a saline calibration tube) sitting in the FOV alongside the leg.
% keep_primary_component tracks the thigh as one connected blob through
% the stack (anchored at the middle slice) and drops everything else.
% See helpers/keep_primary_component.m for the method and its limits.
fprintf('Removing spurious extra tissue blobs (torso, phantom, etc.)...\n');
pd_mask = keep_primary_component(pd_mask);

clear ref_vol n_slices z

if ~exist('output_dir', 'var')
    output_dir = uigetdir(pwd, 'Select output folder to save mask');
end
save(fullfile(output_dir, 'pd_mask_trial_2_otsu_quadrant.mat'), 'pd_mask', '-v7.3');
fprintf('Quadrant-adaptive mask saved to %s\n', output_dir);

