%% 4.1 Denoising with anisotropic smoothing
% Applies edge-preserving anisotropic smoothing to the registered DTI volumes.
% Parameters follow Buck et al. (PLoS ONE, 2015) as used in pre_process.m.
% delta_t is scaled by the actual number of volumes in our data.
% Requires aniso4D_smoothing.m on the MATLAB path (Preprocessing-Functions/).

% trial_1 used noise=5 (sigma=0.05, rho=0.10) and showed tissue-edge
% structure in the Section 4.3 residual map -- a sign of over-smoothing.
% Lower `noise` shrinks sigma/rho together (less edge-blur, less smoothing).
% Change trial_label each time you test a new value so outputs don't
% overwrite each other and Section 4.4 can compare two trials directly.
trial_label = 'trial_2';
noise       = 2.5;
sigma      = noise / 100;
rho        = 2 * sigma;
delta_t    = noise * 3 / size(dti_all_reg, 4);
schemetype = 'implicit_multisplitting';
isnormg    = false;
isfasteig  = true;

% aniso4D_smoothing expects [x_spacing, y_spacing, slice_spacing] (3 elements).
% dti_dims is [in-plane pixel size, slice spacing] — expand to 3-element vector.
% In-plane pixels are isotropic so x and y spacing are the same.
dti_res = [dti_dims(1), dti_dims(1), dti_dims(2)];

fprintf('Running anisotropic smoothing [%s] (noise=%.2f, sigma=%.4f, rho=%.4f, %d volumes, delta_t=%.4f)...\n', ...
    trial_label, noise, sigma, rho, size(dti_all_reg, 4), delta_t);
dti_all_smooth = aniso4D_smoothing(dti_all_reg, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);

diff_max = max(abs(dti_all_smooth(:) - dti_all_reg(:)));
if diff_max < 1e-6
    fprintf('WARNING: Smoothing produced no detectable change (max diff = %.2e).\n', diff_max);
    fprintf('  -> Run: which(''aniso4D_smoothing'') to confirm the function is on the path.\n');
    fprintf('  -> Run: addpath(''tools/MuscleDTI_Toolbox/Preprocessing-Functions'') if needed.\n');
else
    fprintf('Smoothing complete. Max pixel change: %.2f\n', diff_max);
end
clear diff_max

save(fullfile(output_dir, sprintf('dti_smooth_%s_aniso.mat', trial_label)), 'dti_all_smooth', '-v7.3');
fprintf('Saved to %s\n', output_dir);

clear noise sigma rho delta_t schemetype isnormg isfasteig dti_res trial_label

