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
qc_fig_tsnr = figure('Name', 'Denoising Comparison: tSNR', 'Color', 'w', 'Position', [50 50 1400 400]);
subplot(1,4,1); imagesc(tsnr_raw(:,:,z_qc), clim); colormap(gca,'hot'); colorbar;
title(sprintf('Raw (slice %d)', z_qc)); axis image off;
subplot(1,4,2); imagesc(tsnr_a(:,:,z_qc), clim); colormap(gca,'hot'); colorbar;
title(sprintf('Method A — %s', fa)); axis image off;
subplot(1,4,3); imagesc(tsnr_b(:,:,z_qc), clim); colormap(gca,'hot'); colorbar;
title(sprintf('Method B — %s', fb)); axis image off;
subplot(1,4,4); imagesc(tsnr_b(:,:,z_qc) - tsnr_a(:,:,z_qc)); colormap(gca, bwr_colormap()); colorbar;
title('Difference (B − A)'); axis image off;
if exist('output_dir', 'var')
    save_qc_figure(qc_fig_tsnr, output_dir, 'qc_denoising_compare_tsnr');
end

% Signal profile comparison
if exist('pd_mask', 'var')
    [ri, ci] = find(logical(pd_mask(:,:,z_qc)));
    v_row = round(median(ri)); v_col = round(median(ci));
else
    v_row = round(size(dti_all_reg,1)/2); v_col = round(size(dti_all_reg,2)/2);
end

qc_fig_signal = figure('Name', 'Denoising Comparison: Signal Profile', 'Color', 'w', 'Position', [50 100 800 350]);
plot(squeeze(dti_all_reg(v_row,v_col,z_qc,:)), 'b-',  'LineWidth', 1,   'DisplayName', 'Raw'); hold on;
plot(squeeze(smooth_a(v_row,v_col,z_qc,:)),    'r-',  'LineWidth', 1.5, 'DisplayName', fa);
plot(squeeze(smooth_b(v_row,v_col,z_qc,:)),    'g-',  'LineWidth', 1.5, 'DisplayName', fb);
xlabel('Volume index'); ylabel('Signal intensity');
title(sprintf('Signal profile — voxel [%d, %d, slice %d]', v_row, v_col, z_qc));
legend; grid on;
if exist('output_dir', 'var')
    save_qc_figure(qc_fig_signal, output_dir, 'qc_denoising_compare_signal');
end

% Residual comparison
res_a = dti_all_reg(:,:,z_qc,v_qc) - smooth_a(:,:,z_qc,v_qc);
res_b = dti_all_reg(:,:,z_qc,v_qc) - smooth_b(:,:,z_qc,v_qc);
clim_r = [-1 1] * max(max(abs(res_a(:))), max(abs(res_b(:))));
qc_fig_resid = figure('Name', 'Denoising Comparison: Residuals', 'Color', 'w', 'Position', [50 500 800 350]);
subplot(1,2,1); imagesc(res_a, clim_r); colormap(gca, bwr_colormap()); colorbar;
title(sprintf('Residual — %s', fa)); axis image off;
subplot(1,2,2); imagesc(res_b, clim_r); colormap(gca, bwr_colormap()); colorbar;
title(sprintf('Residual — %s', fb)); axis image off;
if exist('output_dir', 'var')
    save_qc_figure(qc_fig_resid, output_dir, 'qc_denoising_compare_residuals');
end

clear tsnr_raw tsnr_a tsnr_b smooth_a smooth_b clim clim_r ...
      res_a res_b v_row v_col z_qc v_qc ri ci fa fb da db start_dir ...
      qc_fig_tsnr qc_fig_signal qc_fig_resid

end

