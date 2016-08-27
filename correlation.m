% correlation matrix
%
sim = all_the_things{1}{2};

wm_ids = [sim.wm_ids(1:4) sim.wm_ids(13:15)];

wm_act = {};
for cond = [1 2]
    A = all_the_things{cond};
    act = A{3};
    onsets = A{4};
    offsets = A{5};
    is_target = A{10};

    pm_onset = onsets(logical(is_target));
    pm_offset = offsets(logical(is_target));
    wm_act{cond} = [];
    for onset = pm_onset
        wm_act_trial = squeeze(act(1, onset:onset+20, wm_ids));
        wm_act{cond} = [wm_act{cond}; wm_act_trial];
    end
end

wm_act_focal = wm_act{1};
wm_act_focal(:,5) = 1e-10; % Monitor tor
wm_act_nonfocal = wm_act{2};
wm_act_nonfocal(:,4) = 1e-10; % Monitor tortoise

n = size(wm_act_focal, 2);
M = corr(wm_act_focal, wm_act_nonfocal);
L = sim.units(wm_ids);

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