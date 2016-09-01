% some knobs
%
save_figures_as_png = true;
save_figures_as_fig = true;

% correlation matrix
%
sim = all_the_things{1}{2};
is_target = all_the_things{1}{10};

% only the "interesting" WM units
wm_ids = [sim.wm_ids(1:4) sim.wm_ids(13:15)];

% pick which trials to analyze
wm_capacity_label = {'WM bias = 4 (med)'};
trial_types = {strcat('OG trials;', {' '}, wm_capacity_label), ...
    strcat('PM trials;', {' '}, wm_capacity_label)};
for PM_trials_or_OG_trials = 0:1 % 1 = PM trials; 0 = OG trials

    trial_type = trial_types{PM_trials_or_OG_trials + 1};
    trials = logical(is_target == PM_trials_or_OG_trials);

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
    block_duration = -1;
    last_block_start = -1;
    for cond = 1:4
        A = all_the_things{cond};
        act = A{3};
        onsets = A{4}(:,trials);
        offsets = A{5}(:,trials);

        wm_act{cond} = [];
        n_subjects = size(act, 1);
        n_trials = size(onsets, 2);
        n_block_trials = n_trials / 4;
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
            % and append them to the WM time series
            wm_act{cond} = [wm_act{cond}; wm_act_trial];
            % mark the boundaries of the first and last blocks
            if i == n_block_trials
                block_duration = size(wm_act{cond}, 1);
            end
            if i + n_block_trials == n_trials
                last_block_start = size(wm_act{cond}, 1) + 1;
            end
        end
    end

    % get the block duration as the min of the first block and the last block
    % durations
    %
    assert(block_duration ~= -1);
    assert(last_block_start ~= -1);
    block_duration = min(size(wm_act{cond}, 1) - last_block_start + 1, block_duration);

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

    % get correlation matrix for condition x vs. condition y
    corr_conds = {};
    corr_conds{1} = {wm_act_focal_low_emph, wm_act_nonfocal_low_emph, ...
        strcat('Focal vs. Nonfocal, Low Emph;', {' '}, trial_type), ...
        'Focal', 'Nonfocal'};
    corr_conds{2} = {wm_act_focal_high_emph, wm_act_nonfocal_high_emph, ...
        strcat('Focal vs. Nonfocal, High Emph;', {' '}, trial_type), ...
        'Focal', 'Nonfocal'};
    corr_conds{3} = {wm_act_focal_low_emph, wm_act_focal_high_emph, ...
        strcat('Low Emph vs. High Emph, Focal;', {' '}, trial_type), ...
        'Low Emph', 'High Emph'};
    corr_conds{4} = {wm_act_nonfocal_low_emph, wm_act_nonfocal_high_emph, ...
        strcat('Low Emph vs. High Emph, Nonfocal;', {' '}, trial_type), ...
        'Low Emph', 'High Emph'};
    corr_conds{5} = {wm_act_focal_low_emph(1:block_duration, :), wm_act_focal_low_emph(end-block_duration+1:end, :), ...
        strcat('Block #1 vs. Block #4, Focal, Low Emph;', {' '}, trial_type), ...
        'Block #1', 'Block #4'};
    corr_conds{6} = {wm_act_focal_high_emph(1:block_duration, :), wm_act_focal_high_emph(end-block_duration+1:end, :), ...
        strcat('Block #1 vs. Block #4, Focal, High Emph;', {' '}, trial_type), ...
        'Block #1', 'Block #4'};
    corr_conds{7} = {wm_act_nonfocal_low_emph(1:block_duration, :), wm_act_nonfocal_low_emph(end-block_duration+1:end, :), ...
        strcat('Block #1 vs. Block #4, Nonfocal, Low Emph;', {' '}, trial_type), ...
        'Block #1', 'Block #4'};
    corr_conds{8} = {wm_act_nonfocal_high_emph(1:block_duration, :), wm_act_nonfocal_high_emph(end-block_duration+1:end, :), ...
        strcat('Block #1 vs. Block #4, Nonfocal, High Emph;', {' '}, trial_type), ...
        'Block #1', 'Block #4'};

    % plot them
    for corr_cond = 1:8
        rows = corr_conds{corr_cond}{1};
        cols = corr_conds{corr_cond}{2};
        fig_title = corr_conds{corr_cond}{3};
        rows_title = corr_conds{corr_cond}{4};
        cols_title = corr_conds{corr_cond}{5};

        n = size(rows, 2);
        M = corr(rows, cols);
        L = sim.units(wm_ids);

        % plot correlation matrix
        figure;
        imagesc(M); % plot the matrix
        set(gca, 'XTick', 1:n); % center x-axis ticks on bins
        set(gca, 'YTick', 1:n); % center y-axis ticks on bins
        set(gca, 'XTickLabel', L); % set x-axis labels
        set(gca, 'YTickLabel', L); % set y-axis labels
        xlabel(cols_title, 'FontSize', 14); % note these are flipped
        ylabel(rows_title, 'FontSize', 14); % note these are flipped
        title(fig_title, 'FontSize', 15); % set title
        colormap('jet'); % set the colorscheme
        colorbar; % enable colorbar
        rotateXLabels(gca, 45);
        
        if save_figures_as_png
            saveas(gcf, strcat('figures/', fig_title{1}, '.png'), 'png');
        end
        if save_figures_as_fig
            saveas(gcf, strcat('figures/', fig_title{1}, '.fig'), 'fig');
        end
    end
    
end % for each trial type (PM / OG)
