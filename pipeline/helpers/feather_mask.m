function soft_mask = feather_mask(mask, sigma)
%FEATHER_MASK  Smooth a binary mask's edges into a gradual taper.
%
% Used before multiplying a mask onto an image for intensity-based
% registration (imregdemons). A hard 0/1 mask boundary creates an
% artificial sharp intensity cliff that demons treats as a strong "edge"
% to align -- if the fixed-side and moving-side masks don't agree exactly
% on where that boundary is (which they generally won't, since they're
% computed independently from different images/contrasts), demons can
% produce large, unstable, oscillatory displacement fields trying to
% reconcile the mismatched hard edges. Feathering removes that false edge
% without requiring the two masks to match exactly.
%
% sigma: Gaussian blur radius in pixels (default 3). Small enough to only
% soften the boundary itself, not blur away real tissue contrast well
% inside the mask.

if nargin < 2
    sigma = 3;
end

soft_mask = imgaussfilt(double(mask), sigma);

end
