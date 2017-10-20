function SBJ08c_HFA_report_actv_cond(SBJ,conditions,an_id,actv_win)
% Load HFA analysis, compute one sample t-test across trials for each time point and electrode,
%   (requires data to be z-scored to baseline), write .txt report
%   of active electrodes and condition differentiating elecs based on stats
% clear all; %close all;

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

stats_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id,'.mat');
load(stats_filename);

% Trim data to stats epoch
cfg_trim = [];
cfg_trim.latency = stat_lim;
hfa{1} = ft_selectdata(cfg_trim,hfa{1});
hfa{2} = ft_selectdata(cfg_trim,hfa{2});
stat   = ft_selectdata(cfg_trim,stat);

% Combine data for stats testing
cfg_comb = [];
cfg_comb.parameter = 'powspctrm';
hfa_all = ft_appendfreq(cfg_comb,hfa{:});

%% Plot Results
% NO! I will plot whatever channels I ran the stats on (why else did I run them?)
% % Select data to plot only ROI channels
% cfgs = [];
% cfgs.channel = {'LPC*'};%SBJ_vars.ch_lab.ROI;
% stat = ft_selectdata(cfgs,stat);
% for an_ix = 1:numel(cond_lab)
%     hfa{an_ix} = ft_selectdata(cfgs,hfa{an_ix});
% end

report_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/' conditions '/' an_id '/'];
if ~exist(report_dir,'dir')
    mkdir(report_dir);
end

% Create a figure for each channel
cond_ch = {};
actv_ch = {};
cond_ch_epochs = {};
actv_ch_epochs = {};
for ch_ix = 1:numel(stat.label)
    figure('Name',stat.label{ch_ix});
    hold on;
    % Compute t-test per time point
    n_tests   = size(hfa_all.powspctrm,4);
    pvals     = NaN([1 n_tests]);
    for t_ix = 1:n_tests
        [~, pvals(t_ix)] = ttest(squeeze(hfa_all.powspctrm(:,ch_ix,1,t_ix)));
    end
    
    % Find epochs with significant task activations
    [~, qvals] = mafdr(pvals);
    actv_mask = qvals<=0.05;
    actv_chunks = fn_find_chunks(actv_mask);
    actv_chunks(actv_mask(actv_chunks(:,1))==0,:) = [];
    actv_chunk_sz = diff(actv_chunks,1,2)+1;
    actv_epochs = actv_chunks(actv_chunk_sz > actv_win,:);
    if ~isempty(actv_epochs)
        actv_ch = {actv_ch{:} stat.label{ch_ix}};
        actv_ch_epochs = {actv_ch_epochs{:}, hfa_all.time(actv_epochs)};
    end
    
%     % Plot distribution of activity
%     edges = [-5:0.2:5];
%     sig_ix = [];
%     for ep_ix = 1:size(actv_epochs,1)
%         sig_ix = [sig_ix actv_epochs(ep_ix,1):actv_epochs(ep_ix,2)];
%     end
%     for t_ix = 1:n_tests
%         if find(sig_ix==t_ix)
%             plot(edges(1:end-1),histcounts(squeeze(hfa_all.powspctrm(:,ch_ix,1,t_ix)),edges),...
%                 'LineWidth',2,'Color','k');
%         else
%             plot(edges(1:end-1),histcounts(squeeze(hfa_all.powspctrm(:,ch_ix,1,t_ix)),edges),...
%                 'LineWidth',1,'Color','c');
%         end
%         title(t_ix);
%         pause(0.1);
%     end
    
    % Find epochs with significant condition differences
    if sum(squeeze(stat.mask(ch_ix,1,:)))>0
        cond_ch = {cond_ch{:} stat.label{ch_ix}};
        mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
        cond_epochs = mask_chunks;
        cond_epochs(squeeze(stat.mask(ch_ix,1,cond_epochs(:,1)))==0,:) = [];
        cond_ch_epochs = {cond_ch_epochs{:} stat.time(cond_epochs)};
        % If stat and hfa aren't on same time axis, adjust sig_chunk indices
        if (size(stat.time,2)~=size(hfa{1}.time,2)) || (sum(stat.time==hfa{1}.time)~=numel(stat.time))
            error('WHy are the stat and hfa on different axes if I selected the data to stat_lim???');
        end
    end
end

%% Print reports of significant channels and their epochs
% Save out list of channels with significant activations from baseline
report_filename = [report_dir 'actv_ch_list.csv'];
report = fopen(report_filename,'w');
fprintf(report,'%s\n',an_id);
for ch_ix = 1:numel(actv_ch)
    for ep_ix = 1:size(actv_ch_epochs{ch_ix},1)
        % Print a row for every significant epoch (multiple rows per channel if needed)
        fprintf(report,'%s,%f,%f\n',actv_ch{ch_ix},actv_ch_epochs{ch_ix}(ep_ix,:));
    end
end
fclose(report);

% Save out list of channels with significant condition differences
report_filename = [report_dir 'cond_ch_list.csv'];
report = fopen(report_filename,'w');
fprintf(report,'%s\n',an_id);
for ch_ix = 1:numel(cond_ch)
    for ep_ix = 1:size(cond_ch_epochs{ch_ix},1)
        % Print a row for every significant epoch (multiple rows per channel if needed)
        fprintf(report,'%s,%f,%f\n',cond_ch{ch_ix},cond_ch_epochs{ch_ix}(ep_ix,:));
    end
end
fclose(report);

end
