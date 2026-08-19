function mask = form_dwi_mask(img)
% FORM_DWI_MASK  Adaptive per-slice mask for a single moving (diffusion-
% weighted) image used in 3_registration/demons.m and demons_parallel.m,
% and in 2_masking/quadrant_otsu.m.
%
% Unlike pd_mask (2_masking/otsu.m or quadrant_otsu.m), which is computed
% once from the anatomical/b0 reference and reused unchanged as the
% FIXED-side mask, this mask is meant to be recomputed fresh for every
% moving slice, from that slice's own contrast. Otsu's threshold is run
% independently per image quadrant (rather than once globally) to tolerate
% signal inhomogeneity across the FOV -- e.g. one leg sitting closer to the
% coil than the other, or general shading -- then the quadrant masks are
% recombined, hole-filled, opened to remove speckle, and dilated once to
% give the tissue edge a small margin. That margin matters here
% specifically: it keeps the true (possibly EPI-distorted) moving tissue
% boundary inside the masked region, so demons still has real signal there
% to align, instead of a hard cutoff that hides the mismatch before
% registration even runs.
%
% Adapted from the quadrant-Otsu masking block in
% bdamon/MuscleDTI_Toolbox Sample-Scripts/pre_process.m (registration
% section), generalized here to non-square matrices and wrapped as a
% reusable function since it's called once per slice/volume rather than
% inline.

[rows, cols] = size(img);
half_r = floor(rows / 2);
half_c = floor(cols / 2);

% quadrant index ranges: [top-left; top-right; bottom-left; bottom-right]
idx_r = [1 half_r; 1 half_r; (half_r + 1) rows; (half_r + 1) rows];
idx_c = [1 half_c; (half_c + 1) cols; 1 half_c; (half_c + 1) cols];

quad_mask = zeros(rows, cols, 4);
for k = 1:4
    r_range = idx_r(k, 1):idx_r(k, 2);
    c_range = idx_c(k, 1):idx_c(k, 2);

    quarter_img = zeros(rows, cols);
    quarter_img(r_range, c_range) = img(r_range, c_range);

    % Guard against an all-zero quadrant (e.g. corner outside the leg)
    quarter_max = max(quarter_img(:));
    if quarter_max > 0
        quarter_img = quarter_img / quarter_max;
    end

    threshold = graythresh(quarter_img); % Otsu's method, this quadrant only
    quarter_mask = zeros(rows, cols);
    quarter_mask(quarter_img > threshold) = 1;
    quad_mask(:, :, k) = quarter_mask;
end

mask = sum(quad_mask, 3) > 0; % recombine the four quadrants
mask = imfill(mask, 'holes');
mask = bwmorph(mask, 'open');   % remove speckle noise
mask = bwmorph(mask, 'dilate'); % one pass -- give the tissue edge margin

end
