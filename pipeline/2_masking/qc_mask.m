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

