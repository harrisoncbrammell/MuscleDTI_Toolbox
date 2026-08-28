%% 3.1 Demons registration method

fprintf('Starting Registration (this may take time)....\n');

% 1. Ensure the mask is a clean matrix
if isstruct(pd_mask) 
    fields = fieldnames(pd_mask);
    pd_mask = pd_mask.(fields{1});
    fprintf("Fixed pd_mask struct. Now it is a matrix.\n");
end

% setup volumes and pre-allocate
[~, ~, slcs, total_vols] = size(dti_all_unreg);
dti_all_reg = zeros(size(dti_all_unreg));

% choose ground truth reference image
if exist('using_anat', 'var') && using_anat && exist('anat_vol', 'var')
    fprintf('Using Anatomical volume as fixed reference for registration.\n');
    fixed_vol = anat_vol;
else
    fprintf('Using average b=0 volume as fixed reference for registration.\n');
    fixed_vol = dti_all_unreg(:, :, :, 1);
    using_anat = false; % enforce false just in case
    dti_all_reg(:, :, :, 1) = dti_all_unreg(:, :, :, 1); % b=0 is perfectly aligned with itself
end

% registration loop slice by slice
fprintf('Starting Slice-by-Slice Demons Registration...\n');

for v = 1:total_vols 
    
    % Skip registering the b=0 image to itself if we aren't using an anatomical reference
    if ~using_anat && v == 1
        continue; 
    end

    fprintf('  Registering Volume %d of %d...\n', v, total_vols);
    
    for z = 1:slcs
        % Grab the slices and current mask
        fixed_slice = fixed_vol(:, :, z);
        moving_slice = dti_all_unreg(:, :, z, v);
        current_mask = pd_mask(:, :, z);

        % Fixed side: mask with the static, pre-computed whole-leg mask.
        % That boundary is trusted -- the anatomical/reference image isn't
        % the one being warped. Feathered (not a hard 0/1 cutoff) because
        % the fixed and moving masks are computed independently and won't
        % agree exactly on where the boundary is -- see feather_mask.m.
        fixed_masked = fixed_slice .* feather_mask(current_mask);

        % Moving side: mask with a FRESH mask built from THIS slice's own
        % contrast, not the shared pd_mask. Reusing pd_mask here would crop
        % the moving image to a boundary computed for a different image,
        % silently assuming the moving tissue edge is already where pd_mask
        % says it should be -- which is exactly what registration is
        % supposed to figure out. Quadrant-adaptive Otsu + dilate gives the
        % moving edge some margin so demons has real signal to align there.
        % See helpers/form_dwi_mask.m; adapted from
        % bdamon/MuscleDTI_Toolbox Sample-Scripts/pre_process.m.
        moving_mask = form_dwi_mask(moving_slice);
        moving_masked = moving_slice .* feather_mask(moving_mask);

        % Safely normalize (Skip registration for this slice if either mask is empty)
        fixed_max = max(fixed_masked(:));
        moving_max = max(moving_masked(:));

        if fixed_max > 0 && moving_max > 0
            fixed_norm = fixed_masked / fixed_max;
            moving_norm = moving_masked / moving_max;
        else
            dti_all_reg(:, :, z, v) = moving_slice; % Copy original and skip
            continue;
        end

        % THE CONTRAST HACK: Apply only when registering b=0 to Anatomical
        if using_anat && v == 1
            moving_norm(moving_norm > 0.5) = 0.5; % Flatten out bright blood vessels
            moving_norm = moving_norm * 2;        % Scale up to match anatomical brightness
        end

        % Calculate displacement field and apply the warp to the UNMASKED original slice
        [disp_field, ~] = imregdemons(moving_norm, fixed_norm, [500 400 200], 'AccumulatedFieldSmoothing', 1.3, 'DisplayWaitbar', false);
        dti_all_reg(:, :, z, v) = imwarp(moving_slice, disp_field);
    end
end

save(fullfile(output_dir, 'dti_registered_trial_1.mat'), 'dti_all_reg', 'pd_mask', 'bvect', '-v7.3');
fprintf('Registration Complete! Saved to %s\n', output_dir);

clear fields rows cols slcs total_vols fixed_vol v z fixed_slice moving_slice current_mask moving_mask fixed_masked moving_masked fixed_max moving_max fixed_norm moving_norm disp_field

