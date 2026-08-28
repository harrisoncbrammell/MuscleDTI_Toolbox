%% 3.4 Cross-Volume Registration QC
%
% Answers: are the nodular artifacts spatially consistent across DWI
% directions, or do they shift per-volume?
%
%   Consistent across volumes → real anatomy (vessels, perimysial septa)
%   Shifting / appearing/disappearing → registration artifact or noise
%
% Outputs two figures:
%   Figure 1 — Montage: a grid of N sampled volumes at one slice, so you
%              can visually scan whether the nodes sit in the same XY
%              position across directions.
%   Figure 2 — Variance map: temporal SD across all registered volumes at
%              that slice. Low SD inside the mask = consistent structure
%              (real). High SD = signal that varies direction-to-direction
%              (could be artifact, noise, or genuine diffusion contrast).

z_inspect = round(size(dti_all_reg, 3) / 2); % middle slice; change as needed
n_show    = 12;   % how many volumes to display in the montage grid

[~, ~, ~, total_vols] = size(dti_all_reg);

% --- Figure 1: Montage across sampled volumes ---

% Sample n_show volumes spread evenly across the direction range,
% skipping vol 1 (b=0 reference -- trivially matched to itself in EPI mode).
vol_idx = round(linspace(2, total_vols, min(n_show, total_vols - 1)));

n_cols = 4;
n_rows = ceil(numel(vol_idx) / n_cols);

fig1 = figure('Name', sprintf('Cross-Volume Montage — Slice %d', z_inspect), ...
               'Color', 'w', 'Position', [50 50 1400 900]);

for k = 1:numel(vol_idx)
    v = vol_idx(k);
    slice_img = dti_all_reg(:, :, z_inspect, v);
    subplot(n_rows, n_cols, k);
    imshow(slice_img, []);
    title(sprintf('Vol %d', v), 'FontSize', 8);
end
sgtitle(sprintf('Registered volumes at slice %d — look for nodes at consistent XY positions', z_inspect), ...
        'FontSize', 10);

% --- Figure 2: Temporal SD map across all volumes ---

% Stack all registered volumes at this slice: [rows x cols x total_vols]
vol_stack = squeeze(dti_all_reg(:, :, z_inspect, :));

% Normalize each volume to [0 1] before computing SD so that global
% intensity differences across b-values don't inflate the variance estimate.
vol_stack_norm = zeros(size(vol_stack));
for v = 1:total_vols
    frame = vol_stack(:, :, v);
    mn = min(frame(:)); mx = max(frame(:));
    if mx > mn
        vol_stack_norm(:, :, v) = (frame - mn) / (mx - mn);
    end
end

sd_map  = std(vol_stack_norm, 0, 3);   % SD across volumes, per pixel
mean_map = mean(vol_stack_norm, 3);     % mean image for context

% Optionally restrict stats to the masked region
if exist('pd_mask', 'var')
    mask_slice = pd_mask(:, :, z_inspect);
else
    mask_slice = ones(size(sd_map));
end

fig2 = figure('Name', sprintf('Temporal SD Map — Slice %d', z_inspect), ...
               'Color', 'w', 'Position', [100 100 1000 420]);

subplot(1, 3, 1);
imshow(mean_map, []);
title('Mean (all vols, normalized)');

subplot(1, 3, 2);
imshow(sd_map, []);
colormap(gca, 'hot'); colorbar;
title('SD across volumes');
xlabel('Low SD = consistent structure | High SD = varies by direction');

subplot(1, 3, 3);
% SD masked to tissue only — makes it easier to see within-muscle variation
sd_masked = sd_map .* mask_slice;
imshow(sd_masked, []);
colormap(gca, 'hot'); colorbar;
title('SD (masked to tissue)');

sgtitle(sprintf('Temporal variance across all %d registered volumes — slice %d', total_vols, z_inspect), ...
        'FontSize', 10);

fprintf('Cross-volume QC done. Slice %d, %d volumes shown in montage.\n', z_inspect, numel(vol_idx));
fprintf('Interpretation:\n');
fprintf('  Montage: if nodes sit at the same XY in every panel -> real anatomy\n');
fprintf('  SD map:  nodes in low-SD regions -> consistent (real); high-SD -> shifts per direction (artifact)\n');

if exist('output_dir', 'var')
    save_qc_figure(fig1, output_dir, sprintf('qc_crossvol_montage_slice%d', z_inspect));
    save_qc_figure(fig2, output_dir, sprintf('qc_crossvol_sdmap_slice%d', z_inspect));
end
