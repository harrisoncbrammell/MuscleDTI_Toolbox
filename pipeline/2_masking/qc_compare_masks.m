%% 2.7 Compare two masks
% Load a second mask and compute Dice similarity — useful for comparing
% masking methods (e.g. Otsu vs define_muscle) or evaluating parameter changes.

start_dir = '';
if exist('output_dir', 'var'), start_dir = output_dir; end

[f1, d1] = uigetfile(fullfile(start_dir, '*.mat'), 'Select FIRST mask');
[f2, d2] = uigetfile(fullfile(start_dir, '*.mat'), 'Select SECOND mask');

if isequal(f1, 0) || isequal(f2, 0)
    fprintf('Comparison cancelled.\n');
else
    m1 = load(fullfile(d1, f1)); if isstruct(m1), fld = fieldnames(m1); m1 = m1.(fld{1}); end
    m2 = load(fullfile(d2, f2)); if isstruct(m2), fld = fieldnames(m2); m2 = m2.(fld{1}); end

    n_slices = size(m1, 3);
    slice_dice = zeros(n_slices, 1);
    for s = 1:n_slices
        a = logical(m1(:,:,s));  b = logical(m2(:,:,s));
        slice_dice(s) = 2*sum(a(:) & b(:)) / (sum(a(:)) + sum(b(:)) + eps);
    end
    overall_dice = 2*sum(m1(:) & m2(:)) / (sum(m1(:)) + sum(m2(:)));

    fprintf('=== Mask Comparison ===\n');
    fprintf('  Mask 1: %s\n', f1);
    fprintf('  Mask 2: %s\n', f2);
    fprintf('  Overall Dice:      %.4f\n', overall_dice);
    fprintf('  Mean per-slice Dice: %.4f\n', mean(slice_dice));
    fprintf('  Min per-slice Dice:  %.4f (slice %d)\n', min(slice_dice), find(slice_dice==min(slice_dice),1));

    qc_fig = figure('Name', 'Mask Comparison: Per-Slice Dice', 'Color', 'w');
    plot(1:n_slices, slice_dice, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.2); hold on;
    yline(overall_dice, 'r--', sprintf('Overall Dice = %.3f', overall_dice), 'LineWidth', 1);
    xlabel('Slice'); ylabel('Dice coefficient');
    title(sprintf('Per-slice Dice: %s vs %s', f1, f2));
    ylim([0 1.05]); grid on;

    if exist('output_dir', 'var')
        save_qc_figure(qc_fig, output_dir, 'qc_mask_compare_dice');
    end

    clear m1 m2 f1 f2 d1 d2 fld n_slices s a b slice_dice overall_dice start_dir qc_fig
end

