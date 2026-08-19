%% 2.5 define_muscle built in masking method

fprintf("Starting manual segmentation using define_muscle...\n");

% Use a temporary reference volume
if exist('using_anat', 'var') && using_anat && exist('anat_vol', 'var')
    fprintf("Using existing anatomical volume for reference.\n");
    ref_vol = anat_vol;
else
    fprintf("No anatomical volume found. Using average b=0 image for anatomical reference.\n");
    ref_vol = dti_all_unreg(:, :, :, 1); 
end

% determine total slices available
total_slices = size(ref_vol, 3);

% prompt user for the slice range they want to segment
prompt = {sprintf('Enter starting slice (1-%d):', total_slices), ...
          sprintf('Enter ending slice (1-%d):', total_slices)};
dlgtitle = 'Slice Range for define_muscle';
dims = [1 50];
definput = {'1', num2str(total_slices)};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    fprintf("Segmentation canceled by user.\n");
else
    slices_to_segment = [str2double(answer{1}), str2double(answer{2})];

    % set up fiber_visualizer options dynamically
    if exist('anat_dims', 'var')
        current_dims = anat_dims;
    else
        current_dims = dti_dims;
    end

    fv_options.anat_dims = current_dims; 
    fv_options.anat_slices = slices_to_segment(1):2:slices_to_segment(2); 
    
    fv_options.plot_mesh = 0;   
    fv_options.plot_mask = 1;   
    fv_options.plot_fibers = 0; 
    
    fv_options.mask_size = [size(ref_vol, 1) size(ref_vol, 2)];
    fv_options.mask_dims = current_dims; % Automatically sets to [FOV SliceThickness]
    fv_options.mask_color = [1 0 0];

    fprintf("Interactive window opening. Follow the define_muscle instructions.\n");
    
    % call the function using our temporary reference volume
    [pd_mask, ~] = define_muscle(ref_vol, slices_to_segment, [], fv_options);
    
    if ~exist('output_dir', 'var')
        output_dir = uigetdir(pwd, 'Select output folder to save mask');
    end
    save(fullfile(output_dir, 'pd_mask_trial_1_define_muscle.mat'), 'pd_mask', '-v7.3');
    fprintf("define_muscle complete! Mask saved to %s\n", output_dir);
end

clear total_slices prompt dlgtitle dims definput answer slices_to_segment fv_options ref_vol current_dims

