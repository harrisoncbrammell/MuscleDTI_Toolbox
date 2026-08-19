%% 5.1 Tensor estimation — weighted least-squares fit (signal2tensor2)
% Fits the Stejskal-Tanner equation per masked voxel using the averaged b=0
% and highest-shell DWI volumes. Follows the pre_process.m toolbox example
% but restructures the loop order (slice outermost) so the inner row/col
% loops run in parallel across cores via parfor.
%
% Inputs required in workspace:
%   dti_all_smooth  (preferred) or dti_all_reg (if denoising was skipped)
%   pd_mask, bvect, bval, highshell_inds, output_dir
%
% Output:
%   tensor_m — [rows × cols × slices × 3 × 3] diffusion tensor field (units: 1/b = mm²/s)

% Use denoised data if available, otherwise fall back to registered data
if exist('dti_all_smooth', 'var')
    dti_input = dti_all_smooth;
    input_label = 'dti_all_smooth';
elseif exist('dti_all_reg', 'var')
    fprintf('WARNING: dti_all_smooth not found. Using undenoised dti_all_reg.\n');
    dti_input = dti_all_reg;
    input_label = 'dti_all_reg';
else
    [f, d] = uigetfile(fullfile(output_dir, '*.mat'), 'Select registered or denoised DTI volume');
    if isequal(f, 0), error('No file selected.'); end
    tmp = load(fullfile(d, f)); flds = fieldnames(tmp);
    dti_input = [];
    for fi = 1:numel(flds)
        if ndims(tmp.(flds{fi})) == 4
            dti_input = tmp.(flds{fi}); break;
        end
    end
    if isempty(dti_input), error('Could not find a 4D array in the selected file.'); end
    input_label = f;
    clear tmp flds f d fi
end

[rows, cols, slices, ~] = size(dti_input);

% Build index vector into dti_input for signal2tensor2.
% Vol 1 is the averaged b=0; the DWI block starts at vol 2.
% highshell_inds indexes within the DWI block (1-based), so +1 shifts
% them to the correct positions in dti_input.
vol_inds = [1; highshell_inds(:) + 1];
n_dirs   = length(highshell_inds);

fprintf('Tensor estimation from: %s\n', input_label);
fprintf('  %d slices, %d DWI directions (b=%.0f s/mm²)\n', slices, n_dirs, bval);
fprintf('  Fitting only within pd_mask (%d masked voxels total)\n', sum(pd_mask(:)));

% Pre-slice into cell arrays — avoids broadcasting the full 4D volume to
% every parfor worker; each worker only receives its own slice.
data_slices = cell(slices, 1);
mask_slices = cell(slices, 1);
for z = 1:slices
    data_slices{z} = dti_input(:, :, z, vol_inds);  % [rows × cols × (1+n_dirs)]
    mask_slices{z} = logical(pd_mask(:, :, z));
end
clear dti_input

% Tensor estimation — parallel over slices
% Each parfor worker fits all voxels within its slice independently.
fprintf('Running tensor estimation (parfor over slices)...\n');
tensor_slices = cell(slices, 1);

parfor z = 1:slices
    slice_data = data_slices{z};    % [rows × cols × (1+n_dirs)]
    slice_mask = mask_slices{z};    % [rows × cols] logical
    t_slice    = zeros(rows, cols, 3, 3);

    for r = 1:rows
        for c = 1:cols
            if slice_mask(r, c)
                signal_v = squeeze(slice_data(r, c, :));
                t_slice(r, c, :, :) = signal2tensor2(signal_v, bvect, bval);
            end
        end
    end

    tensor_slices{z} = t_slice;
    fprintf('  Slice %d/%d done.\n', z, slices);
end

% Reassemble slices into the full 5D tensor volume
tensor_m = zeros(rows, cols, slices, 3, 3);
for z = 1:slices
    tensor_m(:, :, z, :, :) = tensor_slices{z};
end

fprintf('Tensor estimation complete.\n');

save(fullfile(output_dir, 'tensor_trial_1.mat'), 'tensor_m', 'bvect', 'bval', 'pd_mask', '-v7.3');
fprintf('Saved tensor_trial_1.mat to %s\n', output_dir);

clear data_slices mask_slices tensor_slices t_slice slice_data slice_mask ...
      signal_v vol_inds n_dirs rows cols slices z r c input_label

