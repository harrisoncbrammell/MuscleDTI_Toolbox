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
qc_fig = figure('Name', 'Denoising QC', 'Color', 'w', 'Position', [50 50 1100 750]);

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

if exist('output_dir', 'var')
    save_qc_figure(qc_fig, output_dir, 'qc_denoising');
end

clear tsnr_raw tsnr_smooth tsnr_diff clim_tsnr clim_diff residual res_lim ...
      sig_raw sig_smooth v_row v_col z_qc v_qc ri ci smooth_label start_dir ...
      ax1 ax2 ax3 ax4 ax5 ax6 ax7 no_change qc_fig

end

