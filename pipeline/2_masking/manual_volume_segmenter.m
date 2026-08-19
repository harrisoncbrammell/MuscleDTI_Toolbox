%% 2.1 Masking muscle of interest manually with volume segmenter

% Use a temporary reference volume so we don't trick Section 3
if exist('using_anat', 'var') && using_anat && exist('anat_vol', 'var')
    fprintf("Using existing anatomical volume for reference.\n");
    ref_vol = anat_vol;
else
    fprintf("No anatomical volume found. Using average b=0 image for anatomical reference.\n");
    ref_vol = dti_all_unreg(:, :, :, 1); 
end

% Open the app with this data loaded
volumeSegmenter(ref_vol)

clear ref_vol % Keep workspace clean

