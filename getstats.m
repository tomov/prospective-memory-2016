function [OG_RT, OG_RT_SD, OG_Hit, PM_RT, PM_RT_SD, PM_Hit, PM_miss_OG_RT, PM_miss_OG_hit, first_PM_RT] = getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, responses, RTs, act, acc, onsets, offsets, nets, is_target, correct, og_correct, is_inter_task, show_pics, do_print)

n_subjects = size(responses, 1);
n_trials = size(responses, 2);

OG_count = zeros(n_subjects, 1);
PM_count = zeros(n_subjects, 1);

OG_correct_RTs = NaN(n_subjects, n_trials);
PM_hit_RTs = NaN(n_subjects, n_trials); 
false_alarm_RTs = NaN(n_subjects, n_trials); 
OG_wrong_RTs = NaN(n_subjects, n_trials); 
PM_miss_RTs = NaN(n_subjects, n_trials); 
PM_miss_correct_OG_RTs = NaN(n_subjects, n_trials); 
OG_timeout_RTs = NaN(n_subjects, n_trials); 
PM_timeout_RTs = NaN(n_subjects, n_trials); 
first_PM_RT = NaN(n_subjects, 1);

for s=1:n_subjects
    for ord=1:n_trials
        if ~isempty(is_inter_task) && is_inter_task(ord)
            continue
        end
        if strcmp(responses{s, ord}, correct{ord}) == 1
            % right answer
            if is_target(ord) == 0
                % OG correct
                OG_count(s) = OG_count(s) + 1;
                OG_correct_RTs(s, ord) = RTs(s, ord);
            else
                % PM hit
                PM_count(s) = PM_count(s) + 1;
                PM_hit_RTs(s, ord) = RTs(s, ord);
            end
        else
            % wrong answer
            if is_target(ord) == 0
                OG_count(s) = OG_count(s) + 1;
                if strcmp(responses{s, ord}, 'timeout') == 1
                    % timeout
                    OG_timeout_RTs(s, ord) = RTs(s, ord);
                    continue;
                end
                if strcmp(responses{s, ord}, 'PM') == 1
                    % false alarm
                    false_alarm_RTs(s, ord) = RTs(s, ord);
                else
                    % OG wrong
                    OG_wrong_RTs(s, ord) = RTs(s, ord);
                end
            else
                PM_count(s) = PM_count(s) + 1;
                if strcmp(responses{s, ord}, 'timeout') == 1
                    % timeout
                    PM_timeout_RTs(s, ord) = RTs(s, ord);
                    continue;
                end
                % PM miss
                PM_miss_RTs(s, ord) = RTs(s, ord);
                if strcmp(responses{s, ord}, og_correct{ord}) == 1
                    % but still correct OG
                    PM_miss_correct_OG_RTs(s, ord) = RTs(s, ord);
                end
            end
        end
    end
end



%{
RTs
responses
OG_count
PM_count
%}

title_string = 'N/A';
if OG_ONLY
    og_string = 'No PM task';
else
    og_string = 'PM task';
end
if FOCAL
    if EMPHASIS
        title_string = sprintf('focal, high emphasis, %s, %d target(s)', og_string, TARGETS);
    else
        title_string = sprintf('focal, low emphasis, %s, %d target(s)', og_string, TARGETS);
    end
else
    if EMPHASIS
        title_string = sprintf('nonfocal, high emphasis, %s, %d target(s)', og_string, TARGETS);
    else
        title_string = sprintf('nonfocal, low emphasis, %s, %d target(s)', og_string, TARGETS);
    end
end

if ~show_pics && do_print
    for s=1:n_subjects
        fprintf('\n ========================= SUBJECT %d ==========================\n', s);
        if ~show_pics && do_print
            fprintf('\n ----> %s ----\n', title_string);
            fprintf('mean OG correct RTs = %.4f (%.4f)\n', nanmean(OG_correct_RTs(s, :)), nanstd(OG_correct_RTs(s, :)));
            fprintf('mean PM hit RTs = %.4f (%.4f)\n', nanmean(PM_hit_RTs(s, :)), nanstd(PM_hit_RTs(s, :)));
            fprintf('OG accuracy = %.4f%%\n', sum(~isnan(OG_correct_RTs(s, :))) / OG_count(s, :) * 100);
            fprintf('PM hit rate = %.4f%% (%.4f%% were OG correct)\n', sum(~isnan(PM_hit_RTs(s, :))) / PM_count(s, :) * 100, ...
                sum(~isnan(PM_miss_correct_OG_RTs(s, :))) / sum(~isnan(PM_miss_RTs(s, :))) * 100);
        end
    end
end

% return statistics for subject

OG_RT = nanmean(OG_correct_RTs, 2);
% http://en.wikipedia.org/wiki/Standard_error !!!
OG_RT_SD = nanstd(OG_correct_RTs, [], 2) ./ sqrt(sum(~isnan(OG_correct_RTs), 2));
OG_Hit = sum(~isnan(OG_correct_RTs), 2) ./ OG_count * 100;

PM_RT = nanmean(PM_hit_RTs, 2);
for s=1:n_subjects
    for ord=1:n_trials
        if ~isnan(PM_hit_RTs(s, ord))
            first_PM_RT(s) = PM_hit_RTs(s, ord);
            break
        end
    end
    assert(sum(~isnan(PM_hit_RTs(s))) == 0 || ~isnan(first_PM_RT(s)));
end
% http://en.wikipedia.org/wiki/Standard_error !!!
PM_RT_SD = nanstd(PM_hit_RTs, [], 2) ./ sqrt(sum(~isnan(PM_hit_RTs), 2));
PM_Hit = sum(~isnan(PM_hit_RTs), 2) ./ PM_count * 100;

PM_miss_OG_RT = nanmean(PM_miss_correct_OG_RTs, 2);
PM_miss_OG_hit = sum(~isnan(PM_miss_correct_OG_RTs), 2) ./ sum(~isnan(PM_miss_RTs), 2) * 100;

% show figures

act_all = act;
acc_all = acc;
net_all = nets;

if show_pics
    for s=1:n_subjects
        figure;

        act = squeeze(act_all(s,:,:));
        acc = squeeze(acc_all(s,:,:));
        net = squeeze(net_all(s,:,:));

        save('getstats.mat'); % so we can debug more easily; note only last subject is saved

        t_range = 1:min(5000, length(act));
        %t_range = 1:2000;
        y_lim = [sim.MINIMUM_ACTIVATION - 0.1 sim.MAXIMUM_ACTIVATION + 0.1];
        bar_names = {'OG correct', 'PM hit', 'false alarm', 'OG wrong', 'PM miss', 'PM OG' 'OG timeout', 'PM timeout'};
        onset_plot = onsets(onsets < t_range(end))';
        offset_plot = offsets(offsets < t_range(end))';
        % turn off onset plot if necessary
        onset_plot = 0; offset_plot = 0;

        %
        % Experimental stuff
        %

        subplot(4, 3, 1);
        plot(act(1:end, sim.hippo_ids));
        legend(sim.units(sim.hippo_ids));
        title('Hippocampus');
        ylim(y_lim);
        line([onset_plot onset_plot],y_lim,'Color',[0.5 0.5 0.5])
        line([offset_plot offset_plot],y_lim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5])

        subplot(4, 3, 4);
        unit = sim.unit_id('PM Task'); 
        which = sim.weights(:, unit) ~= 0;
        net_inputs = act(:, :) .* repmat(sim.weights(:, unit)', size(act, 1), 1);
        plot(t_range, net_inputs(t_range, which));
        legend(sim.units(which));
        title('Net inputs to PM Task');

        subplot(4, 3, 7);
        unit = sim.unit_id('Monitor tor'); 
        which = sim.weights(:, unit) ~= 0;
        net_inputs = act(:, :) .* repmat(sim.weights(:, unit)', size(act, 1), 1);
        plot(t_range, net_inputs(t_range, which));
        legend(sim.units(which));
        title('Net inputs to Monitor tor');
        %self.net_input(~responded, self.ffwd_ids) = act(:, :) * self.weights(:, self.ffwd_ids) ...
        %    + repmat(self.bias(self.ffwd_ids), sum(~responded), 1);

        subplot(4, 3, 10);
        which = [sim.unit_id('PM Task'), sim.unit_id('hippo 1')];
        plot(t_range, net(t_range, which));
        legend(sim.units(which));
        title('Net input');
        ylim([-10 10]);

        %
        % Stimulus-response pathways
        %
        
        subplot(4, 3, 2);
        plot(t_range, act(t_range, sim.output_ids));
        legend(sim.units(sim.output_ids));
        title('Outputs');
        ylim(y_lim);

        subplot(4, 3, 5);
        plot(t_range, act(t_range, sim.response_ids));
        legend(sim.units(sim.response_ids));
        title('Responses');
        ylim(y_lim);

        subplot(4, 3, 8);
        plot(t_range, act(t_range, sim.perception_ids));
        legend(sim.units(sim.perception_ids));
        title('Feature Perception');
        ylim(y_lim);

        subplot(4, 3, 11);
        plot(t_range, act(t_range, sim.input_ids));
        legend(sim.units(sim.input_ids));
        title('Stimulus Inputs');
        ylim(y_lim);

        %
        % WM submodule & evidence accum
        %

        subplot(4, 3, 3);
        plot(t_range, acc(t_range, :));
        legend(sim.units(sim.output_ids));
        title('Evidence Accumulation');
        %ylim([sim.MINIMUM_ACTIVATION sim.MAXIMUM_ACTIVATION]);

        subplot(4, 3, 6);
        plot(act(1:end, sim.task_ids));
        legend(sim.units(sim.task_ids));
        title('Task Representation');
        ylim(y_lim);
        line([onset_plot onset_plot],y_lim,'Color',[0.5 0.5 0.5])
        line([offset_plot offset_plot],y_lim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5])

        subplot(4, 3, 9);
        plot(act(1:end, sim.attention_ids));
        legend(sim.units(sim.attention_ids));
        title('Feature Attention');
        ylim(y_lim);
        line([onset_plot onset_plot],y_lim,'Color',[0.5 0.5 0.5])
        line([offset_plot offset_plot],y_lim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5])
        
        subplot(4, 3, 12);
        plot(act(1:end, sim.context_ids));
        legend(sim.units(sim.context_ids));
        title('Context');
        ylim(y_lim);
        line([onset_plot onset_plot],y_lim,'Color',[0.5 0.5 0.5])
        line([offset_plot offset_plot],y_lim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5])

        % UG can't move things around with this on top...
%        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%        text(0.5, 1, sprintf('\b%s, subject #%d', title_string, s), 'HorizontalAlignment' ,'center','VerticalAlignment', 'top')

        %
        % -- bar plots
        %{
        figure;

        subplot(1, 2, 1);
        bar([mean(OG_correct_RTs), mean(PM_hit_RTs), ...
            mean(false_alarm_RTs), mean(OG_wrong_RTs), ...
            mean(PM_miss_RTs), mean(PM_miss_correct_OG_RTs), ...
            mean(OG_timeout_RTs), mean(PM_timeout_RTs)]);
        set(gca, 'XTickLabel', bar_names);
        ylim([0 sim.CYCLES_PER_SEC]);
        title('RT (cycles)', 'FontWeight','bold');

        subplot(1, 2, 2);
        bar(100 * [size(OG_correct_RTs, 1) / OG_count, size(PM_hit_RTs, 1) / PM_count, ...
            size(false_alarm_RTs, 1) / OG_count, size(OG_wrong_RTs, 1) / OG_count, ...
            size(PM_miss_RTs, 1) / PM_count, size(PM_miss_correct_OG_RTs, 1) / size(PM_miss_RTs, 1), ...
            size(OG_timeout_RTs, 1) / OG_count, size(PM_timeout_RTs, 1) / PM_count]);
        set(gca, 'XTickLabel', bar_names);
        ylim([0 100]);
        title('Fraction of responses (%)', 'FontWeight','bold');
        ylim([0 100]);
        %}
    end
    
end
