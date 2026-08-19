%% Pipeline Runner
% Runs the full preprocessing pipeline stage by stage. Each stage below
% lists every available method as a `run(...)` call -- exactly one per
% stage should be uncommented (the others are the alternatives you can
% swap in). QC scripts are optional and safe to run in addition to
% whichever method you picked.
%
% Runs in the base workspace (via `run`), same as the old %% cell script --
% every stage still reads/writes the same variables it always did
% (dti_all_reg, pd_mask, dti_all_smooth, etc.), so you can still stop after
% any line and inspect the workspace, or jump in partway if you already
% have a checkpoint .mat file loaded.

pipeline_root = fileparts(mfilename('fullpath'));

addpath(fullfile(pipeline_root, '..', 'Preprocessing-Functions'));
addpath(fullfile(pipeline_root, '..', 'Tractography-Functions'));
addpath(fullfile(pipeline_root, 'helpers'));


%% ===== Stage 1: Ingestion =====

run(fullfile(pipeline_root, '1_ingest', 'load_diffusion.m'));

% Optional -- only needed if you plan to register in anatomical (Dixon)
% space rather than EPI (b=0) space. Sets using_anat = true.
run(fullfile(pipeline_root, '1_ingest', 'load_anatomical.m'));


%% ===== Stage 2: Masking =====
% Choose ONE:

% run(fullfile(pipeline_root, '2_masking', 'manual_volume_segmenter.m'));
% run(fullfile(pipeline_root, '2_masking', 'otsu.m'));
run(fullfile(pipeline_root, '2_masking', 'quadrant_otsu.m'));
% run(fullfile(pipeline_root, '2_masking', 'upload_from_mat.m'));
% run(fullfile(pipeline_root, '2_masking', 'upload_from_nifti.m'));
% run(fullfile(pipeline_root, '2_masking', 'mask_via_define_muscle.m'));

% Optional QC -- safe to run alongside whichever method above:
run(fullfile(pipeline_root, '2_masking', 'qc_mask.m'));
% run(fullfile(pipeline_root, '2_masking', 'qc_compare_masks.m'));  % needs two masks in workspace/on disk


%% ===== Stage 3: Registration =====
% Choose ONE:

% run(fullfile(pipeline_root, '3_registration', 'demons.m'));
run(fullfile(pipeline_root, '3_registration', 'demons_parallel.m'));
% run(fullfile(pipeline_root, '3_registration', 'fsl_eddy.m'));  % stub, not implemented yet

% Optional QC:
run(fullfile(pipeline_root, '3_registration', 'qc_registration.m'));


%% ===== Stage 4: Denoising =====
% Choose ONE:

run(fullfile(pipeline_root, '4_denoising', 'aniso.m'));
% run(fullfile(pipeline_root, '4_denoising', 'tpca.m'));  % stub, not implemented yet

% Optional QC:
run(fullfile(pipeline_root, '4_denoising', 'qc_denoising.m'));
% run(fullfile(pipeline_root, '4_denoising', 'qc_denoising_compare.m'));  % needs two denoised volumes in workspace/on disk


%% ===== Stage 5: Tensor Estimation =====
% Not alternatives -- normally run in order.

run(fullfile(pipeline_root, '5_tensor', 'wls_fit.m'));
run(fullfile(pipeline_root, '5_tensor', 'tensor_maps.m'));
run(fullfile(pipeline_root, '5_tensor', 'residual_map.m'));
