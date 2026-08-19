function show_or_warn(ax, data, clim, cmap, ttl, warn)
% SHOW_OR_WARN  Display an imagesc panel, or a warning text if data is unchanged.
% Used by 4_denoising/qc_denoising.m and qc_denoising_compare.m to handle the
% case where denoising produced no effect.
if warn
    set(ax, 'Color', [0.95 0.95 0.95], 'XColor', 'k', 'YColor', 'k');
    axis(ax, 'off');
    text(ax, 0.5, 0.5, sprintf('No change detected\n(smoothing had no effect)'), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Color', [0.6 0 0], 'FontSize', 10);
else
    imagesc(ax, data, clim); colormap(ax, cmap);
    set(ax, 'Visible', 'off', 'DataAspectRatio', [1 1 1]);
end
title(ax, ttl, 'Color', 'k', 'Visible', 'on', 'FontSize', 9);
end
