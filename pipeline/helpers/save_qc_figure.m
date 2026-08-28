function save_qc_figure(fig, output_dir, name)
%SAVE_QC_FIGURE  Save a QC figure to output_dir as a PNG.
%
% fig: figure handle (e.g. the handle returned by figure(...))
% output_dir: folder to save into
% name: filename without extension, e.g. 'qc_mask_metrics'

if nargin < 3 || isempty(name)
    name = matlab.lang.makeValidName(get(fig, 'Name'));
end

if ~isfolder(output_dir)
    fprintf('WARNING: output_dir does not exist, QC figure "%s" was not saved.\n', name);
    return;
end

out_path = fullfile(output_dir, [name '.png']);
exportgraphics(fig, out_path, 'Resolution', 150);
fprintf('Saved QC figure to %s\n', out_path);

end
