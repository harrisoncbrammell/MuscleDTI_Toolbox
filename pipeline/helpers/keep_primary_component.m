function clean_mask = keep_primary_component(mask_stack, min_area)
%KEEP_PRIMARY_COMPONENT  Remove spurious extra tissue blobs (torso, a QC
% calibration phantom, etc.) from a per-slice mask stack by tracking ONE
% connected component continuously through the volume.
%
% Otsu-style masking only separates tissue from air/background -- it has
% no notion of "muscle" specifically, so anything else radiodense enough
% (torso, a saline phantom taped to the coil, etc.) shows up as its own
% connected component alongside the real thigh. Picking "largest
% component per slice" isn't reliable here since torso's cross-sectional
% area can rival or exceed the thigh's on some slices.
%
% Instead: anchor on the middle slice (assumed clean thigh-only -- take
% its largest component as the starting reference), then walk outward
% slice by slice in both directions. At each step, keep whichever
% candidate component's centroid is closest to the previously-kept
% slice's centroid -- i.e. follow the same physical structure through the
% stack rather than re-deciding "biggest" every slice.
%
% After tracking, a second pass flags slices whose area is a statistical
% outlier (>2SD below the mean of nonzero slices, same test qc_mask.m
% uses) and patches them with pixels both neighboring slices agree are
% tissue -- fixes isolated per-slice Otsu threshold failures without
% touching any of the other slices.
%
% LIMITATION: if the thigh and torso are physically touching with no
% background gap in some slices (e.g. right at the hip), they form ONE
% connected component there and this cannot separate them -- check the
% Section 2.6 QC overlay after running this to confirm those slices look
% right.
%
% mask_stack: logical [rows x cols x slices]
% min_area: components smaller than this (pixels) are treated as noise
%           and ignored (default 200)

if nargin < 2
    min_area = 200;
end

[rows, cols, n_slices] = size(mask_stack);
clean_mask = false(rows, cols, n_slices);

z_anchor = round(n_slices / 2);

anchor_props = regionprops(mask_stack(:,:,z_anchor), 'Area', 'Centroid', 'PixelIdxList');
anchor_props = anchor_props([anchor_props.Area] >= min_area);
if isempty(anchor_props)
    error(['keep_primary_component: no component >= min_area (%d px) found on ' ...
           'anchor slice %d. Pick a different z_anchor or lower min_area.'], min_area, z_anchor);
end
[~, idx] = max([anchor_props.Area]);
anchor_slice = false(rows, cols);
anchor_slice(anchor_props(idx).PixelIdxList) = true;
clean_mask(:,:,z_anchor) = anchor_slice;
anchor_centroid = anchor_props(idx).Centroid;

% walk upward from the anchor to the top slice
prev_centroid = anchor_centroid;
for z = (z_anchor + 1):n_slices
    [clean_mask(:,:,z), prev_centroid] = pick_nearest_component(mask_stack(:,:,z), prev_centroid, min_area, rows, cols);
end

% walk downward from the anchor to slice 1
prev_centroid = anchor_centroid;
for z = (z_anchor - 1):-1:1
    [clean_mask(:,:,z), prev_centroid] = pick_nearest_component(mask_stack(:,:,z), prev_centroid, min_area, rows, cols);
end

% --- Patch per-slice threshold failures (isolated area dips) ---
% Otsu occasionally misses part of the interior on an individual slice
% (local noise/inhomogeneity), which shows up as a random, isolated dip in
% the per-slice area plot rather than a smooth trend. Tuning parameters
% globally to fix a handful of one-off slices risks under-segmenting all
% the other slices that are already fine. Instead, patch only the flagged
% outlier slices by filling in pixels that BOTH neighboring slices agree
% are tissue -- a conservative repair using 3D continuity rather than
% re-thresholding.
areas = squeeze(sum(sum(clean_mask, 1), 2));
nonzero_areas = areas(areas > 0);
if numel(nonzero_areas) > 2
    area_thresh = mean(nonzero_areas) - 2*std(nonzero_areas);
    outlier_z = find(areas > 0 & areas < area_thresh);
    for zi = 1:numel(outlier_z)
        z = outlier_z(zi);
        if z > 1 && z < n_slices
            neighbor_consensus = clean_mask(:,:,z-1) & clean_mask(:,:,z+1);
            clean_mask(:,:,z) = imfill(clean_mask(:,:,z) | neighbor_consensus, 'holes');
        end
    end
end

end


function [slice_mask, centroid] = pick_nearest_component(slice_mask_in, ref_centroid, min_area, rows, cols)
props = regionprops(slice_mask_in, 'Area', 'Centroid', 'PixelIdxList');
props = props([props.Area] >= min_area);

if isempty(props)
    % nothing found on this slice (e.g. past the end of the leg) --
    % return empty and carry the last known centroid forward so tracking
    % can resume if tissue reappears
    slice_mask = false(rows, cols);
    centroid = ref_centroid;
    return;
end

centroids = reshape([props.Centroid], 2, [])';
dists = sqrt(sum((centroids - ref_centroid).^2, 2));
[~, idx] = min(dists);

slice_mask = false(rows, cols);
slice_mask(props(idx).PixelIdxList) = true;
centroid = props(idx).Centroid;
end
