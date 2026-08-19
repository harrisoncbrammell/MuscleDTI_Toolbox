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

