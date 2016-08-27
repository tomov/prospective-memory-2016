% correlation matrix
%
sim = all_the_things{1}{2};

% only the "interesting" WM units
wm_ids = [sim.wm_ids(1:4) sim.wm_ids(13:15)];

% take the shorter duration on each trial (often the focal condition)
% TODO FIXME -- "compress" the nonfocal sequences instead
%
pm_durations = {};
for cond = [1 2]
    A = all_the_things{cond};
    onsets = A{4};
    offsets = A{5};
    is_target = A{10};

    pm_onsets = onsets(logical(is_target));
    pm_offsets = offsets(logical(is_target));
    pm_durations{cond} = pm_offsets - pm_onsets;
end
pm_durations = min(pm_durations{1}, pm_durations{2});

% get the WM activation sequences in the two conditions
wm_act = {};
for cond = [1 2]
    A = all_the_things{cond};
    act = A{3};
    onsets = A{4};
    offsets = A{5};
    is_target = A{10};

    pm_onsets = onsets(logical(is_target));
    pm_offsets = offsets(logical(is_target));
    wm_act{cond} = [];
    for i = 1:size(pm_onsets, 2)
        onset = pm_onsets(i);
        duration = pm_durations(i);
        wm_act_trial = squeeze(act(1, onset:onset+duration-1, wm_ids));
        wm_act{cond} = [wm_act{cond}; wm_act_trial];
    end
end

% a bit of cleanup -- get rid of 0's (they produce NaN's)
wm_act_focal = wm_act{1};
wm_act_focal(:,5) = 1e-10; % Monitor tor
wm_act_nonfocal = wm_act{2};
wm_act_nonfocal(:,4) = 1e-10; % Monitor tortoise

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