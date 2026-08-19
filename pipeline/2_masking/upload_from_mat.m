%% 2.3 Upload mask from mat file
start_dir = '';
if exist('output_dir', 'var'), start_dir = output_dir; end
[mask_file, mask_dir] = uigetfile(fullfile(start_dir, '*.mat'), 'Select the mask .mat file');

if isequal(mask_file, 0)
    fprintf("Mask load cancelled.\n");
else
    pd_mask = load(fullfile(mask_dir, mask_file));
    if isstruct(pd_mask)
        fields = fieldnames(pd_mask);
        pd_mask = pd_mask.(fields{1});
    end
    fprintf("Mask successfully loaded from %s\n", fullfile(mask_dir, mask_file));
end

