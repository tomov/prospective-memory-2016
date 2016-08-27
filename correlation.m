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
    onsets = A{4}(:,trials);
    offsets = A{5}(:,trials);
    durations{cond} = min(offsets - onsets);
    assert(sum(durations{cond} < 0) == 0); % no negative durations allowed
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
    onsets = A{4}(:,trials);
    offsets = A{5}(:,trials);
    
    wm_act{cond} = [];
    n_subjects = size(act, 1);
    n_trials = size(onsets, 2);
    assert(size(onsets, 1) == n_subjects);
    for i = 1:n_trials % for each trial we're interested in
        duration = durations(i);
        wm_act_trial = zeros(n_subjects, duration, length(wm_ids));
        for subj_id = 1:n_subjects % for each subject
            onset = onsets(subj_id, i);
            offset = offsets(subj_id, i);
            % get that subject's WM activations on the trial
            wm_act_trial(subj_id,:,:) = squeeze(act(subj_id, onset:onset + duration - 1, wm_ids));
            %wm_act_trial(subj_id,:,:) = squeeze(act(subj_id, offset - duration + 1:offset, wm_ids));
        end
        % then average the trial WM activations across all subjects
        wm_act_trial = squeeze(mean(wm_act_trial, 1));
        % and append them to the time series
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