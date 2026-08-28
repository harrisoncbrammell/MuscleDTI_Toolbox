%% 4.2 Denoising with TPCA (Threshold PCA)
% Applies TPCA denoising to the registered DTI volumes as an alternative
% to anisotropic smoothing. Removes thermal noise by exploiting the low-rank
% structure of the diffusion signal across gradient directions.
% Requires TPCA_denoising.m on the MATLAB path (Preprocessing-Functions/).
%
% Reference: Henriques et al., Imaging Neuroscience, 2023.
%
% Inputs expected in workspace:
%   dti_all_reg  — [x, y, z, N] registered DTI volumes
%   pd_mask      — [x, y, z] binary tissue mask
%   output_dir   — path string for saving outputs

trial_label = 'trial_1';

% Kernel size for the sliding PCA patch window.
% Default is [5 5 5] but our slice thickness (4mm) is twice the in-plane
% resolution (2mm), so a [5 5 3] kernel is more physically isotropic.
kernel = [5 5 3];


%% Step 1 — Identify background voxels
mean_signal = mean(dti_all_reg, 4);
background_mask = ~pd_mask & (mean_signal < 25);
background_mask = ~imfill(~background_mask, 'holes');
for z = 1:size(background_mask, 3)
    background_mask(:,:,z) = bwareafilt(background_mask(:,:,z), 1); %keep single largest white (background) region per slice
    background_mask(:,:,z) = ~bwareaopen(~background_mask(:,:,z), 50); %fills small holes in mask of up to 50 pixels in area
    background_mask(:,:,z) = imerode(background_mask(:,:,z), strel('disk', 15)); %erode the mask to make sure it doesnt include tissue boarders
end

fprintf('Background voxels: %d\n', sum(background_mask(:)));
z = round(size(background_mask,3)/2);
qc_fig = figure;
subplot(1,2,1); imshow(dti_all_reg(:,:,z,1) .* double(background_mask(:,:,z)), []); title('Masked Image');
subplot(1,2,2); imshow(background_mask(:,:,z)); title('Mask');
if exist('output_dir', 'var'), save_qc_figure(qc_fig, output_dir, 'qc_tpca_background_mask'); end

%% Step 2 — Compute noise variance per background voxel
background_var_map = var(dti_all_reg, 0, 4) .* double(background_mask);
background_var_map(~background_mask) = median(background_var_map(background_mask > 0), 'all');
disp(median(background_var_map(background_mask > 0), 'all'));
%mabye try not zeroing out the tissue?

%% Step 3 — Smooth the noise variance map to fill in masked regions
full_var_map = imgaussfilt3(background_var_map, 50);
disp(median(full_var_map(pd_mask > 0), 'all'));
disp(max(background_var_map(:)));
disp(sum(background_mask(:)));

%% Step 4 — Call TPCA_denoising
[dti_all_smooth, n_comps] = TPCA_denoising(dti_all_reg, pd_mask, kernel, full_var_map);

%% Step 5 — Save
save(fullfile(output_dir, sprintf('dti_smooth_%s_tpca.mat', trial_label)), 'dti_all_smooth', '-v7.3');
fprintf('Saved to %s\n', output_dir);



% Cleanup
clear trial_label kernel background_mask noise_var n_comps
