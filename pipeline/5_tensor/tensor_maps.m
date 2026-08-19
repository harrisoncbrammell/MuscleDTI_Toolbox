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

