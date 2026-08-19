%% 3.3 Registration Quality Control (Visual Inspection)

% 1. Settings: Choose which volume and slice to inspect
% (Volume 2 is usually the first diffusion direction)
v_inspect = 1; % b=0 matches Dixon water contrast best — most diagnostic; change to 2+ for a DWI direction
z_inspect = round(size(dti_all_reg, 3) / 2); % Middle slice

% 2. Prepare the images
% Grab the fixed reference used in Section 3.1
if exist('using_anat', 'var') && using_anat
    fixed_img = anat_vol(:,:,z_inspect);
    ref_label = 'Anatomical';
else
    fixed_img = dti_all_unreg(:,:,z_inspect, 1);
    ref_label = 'Average b=0';
end

moving_unreg = dti_all_unreg(:,:,z_inspect, v_inspect);
moving_reg = dti_all_reg(:,:,z_inspect, v_inspect);

% 3. Create the QC Figure
figure('Name', 'Registration Quality Control', 'NumberTitle', 'off', 'Color', 'w');

% Subplot 1: Checkerboard (Alignment Test)
subplot(2,2,1);
imshowpair(fixed_img, moving_reg, 'checkerboard');
title(['Checkerboard: ', ref_label, ' vs Reg Vol ', num2str(v_inspect)]);
xlabel('Discontinuities in lines = Bad Alignment');

% Subplot 2: Difference Map (Error Heatmap)
% We normalize both to [0 1] before subtracting to see pure alignment errors
f_norm = fixed_img / max(fixed_img(:));
m_norm = moving_reg / max(moving_reg(:));
diff_map = abs(f_norm - m_norm);

subplot(2,2,2);
imshow(diff_map, []);
colormap(gca, 'hot'); 
colorbar;
title('Difference Map (Normalized)');
xlabel('Bright outlines = Registration Error');

% Subplot 3: Overlay (Green/Magenta)
% Fixed is Green, Moving is Magenta. Grey areas are perfectly aligned.
subplot(2,2,3);
imshowpair(fixed_img, moving_reg, 'falsecolor');
title('Overlay: Green(Fixed) / Mag(Reg)');
xlabel('Color fringes indicate mismatch');

% Subplot 4: Before vs After (Side-by-Side)
subplot(2,2,4);
imshowpair(moving_unreg, moving_reg, 'montage');
title('Before (Left) vs After (Right)');
xlabel('Compare distortion/warping');

fprintf('QC Figure generated for Slice %d, Volume %d.\n', z_inspect, v_inspect);

% 4. Optional: Animation (The Flicker Test)
fprintf('Starting Flicker Test in a new window...\n');
figure('Name', 'Flicker Test: Toggle between Fixed and Registered', 'Color', 'k');
for i = 1:6
    imshow(fixed_img, []); title(['FIXED (', ref_label, ')']); pause(0.4);
    imshow(moving_reg, []); title(['REGISTERED (Vol ', num2str(v_inspect), ')']); pause(0.4);
end

