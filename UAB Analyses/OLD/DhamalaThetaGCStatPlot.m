thetacon = squeeze(mean(con(frex >= 4 & frex <= 8,:,:)));
thetainc = squeeze(mean(inc(frex >= 4 & frex <= 8,:,:)));


[hl,~] = boundedline(time,mean(thetainc-thetacon,2),std(thetainc-thetacon,[],2)./sqrt(size(thetacon,2)),'r');
%% Early
stat_con = squeeze(mean(mean(con(frex >= 4 & frex <= 8,time >= 0.25 & time <= 0.526,:)),2));
stat_inc = squeeze(mean(mean(inc(frex >= 4 & frex <= 8,time >= 0.25 & time <= 0.526,:)),2));

% Number of data points
n = length(stat_con);

% X positions for the two conditions
x_con = ones(1, n) + randn(1, n) * 0.025; % Adding jitter to x = 1
x_inc = ones(1, n) * 2 + randn(1, n) * 0.025; % Adding jitter to x = 2

% Define colors
color_con = [0.75, 0.75, 1]; 
color_inc = [1, 0.75, 0.75]; 

% Create figure
figure;
hold on;

% Plot black lines connecting the pairs
for i = 1:n
    plot([x_con(i), x_inc(i)], [stat_con(i), stat_inc(i)], 'k-', 'LineWidth', 0.1);
end

% Plot the circles for cI condition (x = 1)
scatter(x_con, stat_con, 100, color_con, 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1);

% Plot the circles for iI condition (x = 2)
scatter(x_inc, stat_inc, 100, color_inc, 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1);

% Customize plot
set(gca, 'XTick', [1 2], 'XTickLabel', {'Con', 'Inc'});
xlim([0 3]);
ylabel('Net dlPFC -> dACC Theta GC');
title('Early Trial Theta GC');
ylim([-0.06 0.06])
yticks([-0.04 0 0.04])
hold off;

%% Late trial
% Number of data points
stat_con = squeeze(mean(mean(con(frex >= 4 & frex <= 8,time >= 0.65 & time <= 0.874,:)),2));
stat_inc = squeeze(mean(mean(inc(frex >= 4 & frex <= 8,time >= 0.65 & time <= 0.874,:)),2));
n = length(stat_con);

% X positions for the two conditions
x_con = ones(1, n) + randn(1, n) * 0.025; % Adding jitter to x = 1
x_inc = ones(1, n) * 2 + randn(1, n) * 0.025; % Adding jitter to x = 2

% Define colors
color_con = [0.75, 0.75, 1]; 
color_inc = [1, 0.75, 0.75]; 

% Create figure
figure;
hold on;

% Plot black lines connecting the pairs
for i = 1:n
    plot([x_con(i), x_inc(i)], [stat_con(i), stat_inc(i)], 'k-', 'LineWidth', 0.1);
end

% Plot the circles for Con condition (x = 1)
scatter(x_con, stat_con, 100, color_con, 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1);

% Plot the circles for Inc condition (x = 2)
scatter(x_inc, stat_inc, 100, color_inc, 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1);


y_max = max([stat_con; stat_inc]); 
y_line = y_max + 0.005; 
plot([1, 2], [y_line, y_line], 'k-', 'LineWidth', 1); 


x_mid = (1 + 2) / 2; 
plot(x_mid, y_line + 0.0015, 'k*', 'MarkerSize', 4.5, 'LineWidth', 1.0); 

% Customize plot
set(gca, 'XTick', [1 2], 'XTickLabel', {'Con', 'Inc'});
xlim([0 3]);
ylabel('Net dlPFC -> dACC Theta GC');
title('Late Trial Theta GC');
ylim([-0.045 0.045]) 
yticks([-0.04 0 0.04])
hold off;
