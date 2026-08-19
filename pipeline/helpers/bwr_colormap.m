function cmap = bwr_colormap()
% BWR_COLORMAP  Blue-white-red diverging colormap (replacement for 'rdbu'
% which is not a built-in MATLAB colormap). Blue = negative, white = zero,
% red = positive. Used by 4_denoising/qc_denoising.m and qc_denoising_compare.m.
n = 64;
cmap = [linspace(0,1,n)' linspace(0,1,n)' ones(n,1); ...
        ones(n,1) linspace(1,0,n)' linspace(1,0,n)'];
end
