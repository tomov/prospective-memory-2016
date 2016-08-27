% correlation matrix
%
sim = all_the_things{1}{2};
is_target = all_the_things{1}{10};

% only the "interesting" WM units
wm_ids = [sim.wm_ids(1:4) sim.wm_ids(13:15)];

% which trials to analyze
trials = logical(~is_target);

% take the shorter duration on each trial (often the focal condition)
% TODO FIXME -- "compress" the nonfocal sequences instead
%
durations = {};
for cond = 1:4
    A = all_the_things{cond};
    onsets = A{4}(trials);
    offsets = A{5}(trials);
    durations{cond} = offsets - onsets;
end
x = durations{1};
for cond = 3
    x = min(x, durations{cond});
end
durations = x;

% get the WM activation sequences in the two conditions
wm_act = {};
for cond = 1:4
    A = all_the_things{cond};
    act = A{3};
    onsets = A{4}(trials);
    offsets = A{5}(trials);
    
    wm_act{cond} = [];
    for i = 1:size(onsets, 2)
        onset = onsets(i);
        duration = durations(i);
        wm_act_trial = squeeze(act(1, onset:onset+duration-1, wm_ids));
        wm_act{cond} = [wm_act{cond}; wm_act_trial];
    end
end

% get the wm activations in the different conditions
% order is determined in EM2005.m -- make sure to use for, not parfor for
% looping over the conditions there
assert(all_the_things{1}{1}(2) == 1 && all_the_things{1}{1}(3) == 0);
assert(all_the_things{2}{1}(2) == 1 && all_the_things{2}{1}(3) == 1);
assert(all_the_things{3}{1}(2) == 0 && all_the_things{3}{1}(3) == 0);
assert(all_the_things{4}{1}(2) == 0 && all_the_things{4}{1}(3) == 1);
wm_act_focal_low_emph = wm_act{1};
wm_act_focal_high_emph = wm_act{2};
wm_act_nonfocal_low_emph = wm_act{3};
wm_act_nonfocal_high_emph = wm_act{4};

% get correlation matrix
n = size(wm_act_focal_low_emph, 2);
M = corr(wm_act_focal_low_emph, wm_act_nonfocal_low_emph);
L = sim.units(wm_ids);

% plot correlation matrix
figure;
imagesc(M); % plot the matrix
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
set(gca, 'XTickLabel', L); % set x-axis labels
set(gca, 'YTickLabel', L); % set y-axis labels
xlabel('Nonfocal', 'FontSize', 14);
ylabel('Focal', 'FontSize', 14);
title('Focal vs. Nonfocal', 'FontSize', 15); % set title
colormap('jet'); % set the colorscheme
colorbar; % enable colorbar