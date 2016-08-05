% focal
A = all_the_things{1};

sim = A{2};
act = A{3};
onsets = A{4};
offsets = A{5};
is_target = A{10};


focal_pm_onset = onsets(logical(is_target));
focal_pm_offset = offsets(logical(is_target));
focal_pm_act = squeeze(act(1, focal_pm_onset(1):focal_pm_offset(1), :));
focal_pm_nodes = mean(focal_pm_act);


figure;
barweb(focal_pm_nodes, std(focal_pm_act)/sqrt(size(focal_pm_act, 1)), 1, {}, ...
    'Focal, Low Emph, PM trial', 'Node', 'Average Activation');
h = legend(sim.units);
set(h, 'FontSize', 10);


% nonfocal
A = all_the_things{2};

sim = A{2};
act = A{3};
onsets = A{4};
offsets = A{5};
is_target = A{10};

nonfocal_pm_onset = onsets(logical(is_target));
nonfocal_pm_offset = offsets(logical(is_target));
nonfocal_pm_act = squeeze(act(1, nonfocal_pm_onset(1):nonfocal_pm_offset(1), :));
nonfocal_pm_nodes = mean(nonfocal_pm_act);


d = pdist2(focal_pm_nodes', nonfocal_pm_nodes');

figure;
imshow(1 - d);


figure;
barweb(nonfocal_pm_nodes, std(nonfocal_pm_act)/sqrt(size(nonfocal_pm_act, 1)), 1, {}, ...
    'Nonfocal, Low Emph, PM trial', 'Node', 'Average Activation');
h = legend(sim.units);
set(h, 'FontSize', 10);


[X,Y] = meshgrid(-8:.5:8);
R = sqrt(X.^2 + Y.^2) + eps;
Z = sin(R)./R;

% surface in 3D
figure;
surf(Z,'EdgeColor','None');

% 2D map using view
figure;
surf(Z,'EdgeColor','None');
view(2);    

% using imagesc to view just Z
figure;
imagesc(Z); 
colormap jet;

% smooth it out
figure;
surf(X, Y, Z,'EdgeColor', 'None', 'facecolor', 'interp');
view(2);
axis equal; 
axis off;

% colorbar

colorbarexample
colorbar(placement)example
colorbar(Name,Value)example
colorbar(placement,Name,Value)
colorbar(ax,___)
colorbar('peer',ax___)
h = colorbar(___)example
colorbar('off')example
colorbar(h,'off')
colorbar(ax,'off')