
x = [-1:0.2:5];
x = 1;
cycs = 250;

%{
I = 0;
li = -2; % lateral inhibition
se = -2; % self-excitation
ve = 0; % vertical excitation
ci = 0; % cross-inhibition
step = 0.2;
%}

I = 2;
li = -2; % lateral inhibition
se = -2 + 0.001; % self-excitation
ve = -1; % vertical excitation
ci = -1; % cross-inhibition
step = 0.1;


W = [
    se li ve ci 0;
    li se ci ve 0;
    ve ci se li 0;
    ci ve li se 0;
    0 0 0 0 0;
    ];
b_high = [I I 0 0 0];
b_low = b_high / 10;
b = b_high;
init_a = [0.99 0.6 0.4 0.3 0];

z = zeros(size(x, 2));

for i = 1:size(x, 2)
    for j = 1:size(x, 2)
        a1 = x(i);
        a2 = x(j);
 
        a = [a1 a2 a1 a2];
        if size(x, 2) == 1
            a = init_a;
        end
        net = [0 0 0 0];
        avg = [0 0 0 0];
        act = zeros(cycs+1, 5);
                                                %act(1,:) = 1 ./ (1 + exp(-a));
                                                act(1,:) = a;
        for cyc = 1:cycs
            net = a * W + b;
            
            if cyc > 90 && cyc < 130
                b(2) = 4;
            else
                if cyc == 130 || cyc == 131
                    a(1:2) = init_a(1:2);
                end
                b(2) = I;
            end
            
            
            %b = b - b * 0.00001;
            
            %{
            if cyc > 50 && cyc < 70
                b = b_low;
            else
                b = b_high;
            end
            %}

            
            %if cyc > 20 && cyc < 40
            %    b(5) = 15;
            %else
            %    b(5) = I;
            %end
            
                                                da = net;
                                                %da = -a + quadsquare(b + a * W);
            a = a + step * da;
            a = min(a, 1);
            a = max(a, 0);
                                                %act(cyc+1,:) = 1 ./ (1 + exp(-a));
                                                act(cyc+1,:) = a;
        end
        z(i, j) = norm(act(end,:) - act(1,:));
        
        if size(x, 2) == 1
            figure;
            
            plot(act(:,1:2), 'LineWidth', 2);
            ylim([0 1]);
            xlim([0 210]);
            legend('OG Task', 'PM Task', 'OG Features', 'PM Features', 'Load');
            ylabel('activation');
            xlabel('cycle');
            line([10 10],ylim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5]);
            line([50 50],ylim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5]);
            line([90 90],ylim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5]);
            
            line([130 130],ylim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5]);
            line([140 140],ylim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5]);
            line([180 180],ylim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5]);
            
            text(12, 0.1, 'trial #1');
            text(52, 0.1, 'trial #2');
            
            text(92, 0.1, 'trial #3');
            text(132, 0.15, 'end trial');            
            text(142, 0.1, 'trial #4');
            text(182, 0.1, 'trial #5');
            
            %title('Ongoing Trials, 5 WM units');            
            %title('Task Switch, 5 WM units');
            
            figureHandle = gcf;
            set(findall(figureHandle,'type','text'),'fontSize',14);
        end
    end
end

if size(x, 2) > 1
    figure;
    surf(x, x, z);
end