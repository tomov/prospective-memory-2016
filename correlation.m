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
for cond = [1 2]
    A = all_the_things{cond};
    onsets = A{4}(trials);
    offsets = A{5}(trials);
    durations{cond} = offsets - onsets;
end
durations = min(durations{1}, durations{2});

% get the WM activation sequences in the two conditions
wm_act = {};
for cond = [1 2]
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

% a bit of cleanup -- get rid of 0's (they produce NaN's)
wm_act_focal = wm_act{1};
%wm_act_focal(:,5) = 1e-10; % Monitor tor
wm_act_nonfocal = wm_act{2};
%wm_act_nonfocal(:,4) = 1e-10; % Monitor tortoise

% get correlation matrix
n = size(wm_act_focal, 2);
M = corr(wm_act_focal, wm_act_nonfocal);
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