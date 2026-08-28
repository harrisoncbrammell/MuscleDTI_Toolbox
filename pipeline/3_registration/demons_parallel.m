%% 3.1.1 Demons registration method cpu accelorated
fprintf('Starting Registration (this may take time)....\n');

% ensure the mask is a clean matrix
if isstruct(pd_mask) 
    fields = fieldnames(pd_mask);
    pd_mask = pd_mask.(fields{1});
    fprintf("Fixed pd_mask struct. Now it is a matrix.\n");
end

% setup volumes and pre-allocate
[rows, cols, slcs, total_vols] = size(dti_all_unreg);
dti_all_reg = zeros(size(dti_all_unreg));

% pick ground truth reference image
if exist('using_anat', 'var') && using_anat && exist('anat_vol', 'var')
    fprintf('Using Anatomical volume as fixed reference for registration.\n');
    fixed_vol = anat_vol;
else
    fprintf('Using average b=0 volume as fixed reference for registration.\n');
    fixed_vol = dti_all_unreg(:, :, :, 1);
    using_anat = false; % enforce false just in case
    dti_all_reg(:, :, :, 1) = dti_all_unreg(:, :, :, 1); % b=0 is perfectly aligned with itself
end

% Pre-slice fixed_vol and pd_mask into cell arrays so parfor workers receive
% only the slices they need rather than broadcasting the entire 3D volume.
fixed_slices = cell(slcs, 1);
mask_slices  = cell(slcs, 1);
for z = 1:slcs
    fixed_slices{z} = fixed_vol(:, :, z);
    mask_slices{z}  = pd_mask(:, :, z);
end

% The Parallel Registration Loop
fprintf('Starting Parallel Slice-by-Slice Demons Registration...\n');
fprintf('Booting up parallel pool (this takes a few seconds the first time)...\n');

parfor v = 1:total_vols 
    
    % Skip registering the b=0 image to itself if we aren't using an anatomical reference
    if ~using_anat && v == 1
        continue; 
    end

    fprintf('  Core processing Volume %d of %d...\n', v, total_vols);
    
    % Create a temporary 3D volume for this specific worker core to build
    temp_vol_reg = zeros(rows, cols, slcs);
    
    for z = 1:slcs
        % Use pre-sliced cell arrays to avoid broadcasting full 3D volumes
        fixed_slice  = fixed_slices{z};
        moving_slice = dti_all_unreg(:, :, z, v);
        current_mask = mask_slices{z};

        % Fixed side: mask with the static, pre-computed whole-leg mask.
        % That boundary is trusted -- the anatomical/reference image isn't
        % the one being warped. Feathered (not a hard 0/1 cutoff) because
        % the fixed and moving masks are computed independently and won't
        % agree exactly on where the boundary is -- see feather_mask.m.
        fixed_masked = fixed_slice .* feather_mask(current_mask);

        % Moving side: mask with a FRESH mask built from THIS slice's own
        % contrast, not the shared pd_mask. See the matching comment in
        % demons.m and helpers/form_dwi_mask.m for the rationale.
        % form_dwi_mask and feather_mask are both self-contained (stateless
        % besides their own inputs), so they're safe to call inside parfor.
        moving_mask = form_dwi_mask(moving_slice);
        moving_masked = moving_slice .* feather_mask(moving_mask);

        % Safely normalize (Skip registration for this slice if either mask is empty)
        fixed_max = max(fixed_masked(:));
        moving_max = max(moving_masked(:));

        if fixed_max > 0 && moving_max > 0
            fixed_norm = fixed_masked / fixed_max;
            moving_norm = moving_masked / moving_max;
        else
            temp_vol_reg(:, :, z) = moving_slice; % Copy original and skip
            continue;
        end

        % THE CONTRAST HACK: Apply only when registering b=0 to Anatomical
        if using_anat && v == 1
            moving_norm(moving_norm > 0.5) = 0.5; % Flatten out bright blood vessels
            moving_norm = moving_norm * 2;        % Scale up to match anatomical brightness
        end

        % Calculate displacement field and apply the warp to the UNMASKED original slice
        [disp_field, ~] = imregdemons(moving_norm, fixed_norm, [500 400 200], 'AccumulatedFieldSmoothing', 2.5, 'DisplayWaitbar', false);
        temp_vol_reg(:, :, z) = imwarp(moving_slice, disp_field);
    end

    % Drop the finished volume into the master array
    dti_all_reg(:, :, :, v) = temp_vol_reg;
end

save(fullfile(output_dir, 'dti_registered_trial_1.mat'), 'dti_all_reg', 'pd_mask', 'bvect', '-v7.3');
fprintf('Registration Complete! Saved to %s\n', output_dir);

clear fields rows cols slcs total_vols fixed_vol v z fixed_slice moving_slice current_mask moving_mask fixed_masked moving_masked fixed_max moving_max fixed_norm moving_norm disp_field temp_vol_reg
