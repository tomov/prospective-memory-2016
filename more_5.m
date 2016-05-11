%{
 subjects[subj_id, :] ==
 samples variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - EMPHASIS
 statistics order = 
    4 - OG_RT, 
    5 - OG_Hit, 
    6 - PM_RT, 
    7 - PM_Hit,
    8 - PM_miss_OG_hit
    9 - third task target RT's
    10 - third task nontarget RT's
    11 - third task hit rate
 (see EM2005)
%}



Ms = [];
SEMs = [];



GroupA = data(:, 9);
%GroupA = GroupA(data(:,1) == 0 & data(:, 3) == 0);

m = mean(GroupA) * RT_slope + RT_intercept
s = std(GroupA) * RT_slope / sqrt(length(GroupA))

Ms = [Ms; m];
SEMs = [SEMs; s];

GroupB = data(:, 10);
%GroupB = GroupB(data(:,1) == 0 & data(:, 3) == 1);

m = mean(GroupB) * RT_slope + RT_intercept
s = std(GroupB) * RT_slope / sqrt(length(GroupB))

Ms = [Ms; m];
SEMs = [SEMs; s];

all = [GroupA; GroupB];
groups = [repmat({'1'}, length(GroupA), 1); repmat({'4'}, length(GroupB), 1)];

[p, table] = anova1(all, groups, 'off')





figure;
barweb(Ms, SEMs, 1, {}, ...
    'Simulation Data', 'Third Task Trial Type');
h = legend({'Target', 'Nontarget'});
set(h, 'FontSize', 15);
ylabel('Third Task RT (ms)');
ylim([900, 1000]);
        

        


%{
RT_slope = 12.5;
RT_intercept = 63;

GroupA = data(:, 5);
%GroupA = GroupA(data(:,1) == 0 & data(:, 3) == 0);

%mean(GroupA) * RT_slope + RT_intercept
%std(GroupA) * RT_slope / sqrt(length(GroupA))
mean(GroupA)
std(GroupA) / sqrt(length(GroupA))

GroupB = data(:, 11);
%GroupB = GroupB(data(:,1) == 0 & data(:, 3) == 1);

%mean(GroupB) * RT_slope + RT_intercept
%std(GroupB) * RT_slope / sqrt(length(GroupB))
mean(GroupB)
std(GroupB) / sqrt(length(GroupB))


all = [GroupA; GroupB];
groups = [repmat({'1'}, length(GroupA), 1); repmat({'4'}, length(GroupB), 1)];

[p, table] = anova1(all, groups, 'off')
%}