function SBJ08b_HFA_plot_SR_actv(SBJ,an_id_s,an_id_r,actv_win,plt_id,save_fig,fig_vis,fig_ftype)
% Plots both stimulus- and response-locked HFA computed in SBJ08ab_HFA_actv
% clear all; %close all;

if ischar(save_fig); save_fig = str2num(save_fig); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set up paths
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

evnt_lab = {'S', 'R'};
% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
srate = trial_info.sample_rate;

hfa_fname1 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_s,'.mat');
hfa_fname2 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_r,'.mat');
tmp = load(hfa_fname1,'hfa'); hfa{1} = tmp.hfa;
tmp = load(hfa_fname2,'hfa'); hfa{2} = tmp.hfa;
tmp_an1 = load(hfa_fname1,'an'); tmp_an2 = load(hfa_fname2,'an');
if any(~[strcmp(tmp_an1.evnt_lab,'S') strcmp(tmp_an2.evnt_lab,'R')]);
    error('HFA analyses not ordered S then R!');
end
stat_fname1 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_s,'_actv_mn',actv_win,'.mat');
stat_fname2 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_r,'_actv_mn',actv_win,'.mat');
tmp = load(stat_fname1,'actv_ch'); actv_ch{1} = tmp.actv_ch;
tmp = load(stat_fname2,'actv_ch'); actv_ch{2} = tmp.actv_ch;
tmp = load(stat_fname1,'actv_ch_epochs'); actv_ch_epochs{1} = tmp.actv_ch_epochs;
tmp = load(stat_fname2,'actv_ch_epochs'); actv_ch_epochs{2} = tmp.actv_ch_epochs;
tmp_st1 = load(stats_fname1,'st'); tmp_st2 = load(stats_fname2,'st');
if any([~strcmp(tmp_an1.evnt_lab,tmp_st1.evnt_lab) ~strcmp(tmp_an2.evnt_lab,tmp_st2.evnt_lab)]);
    error('an and st evnt_lab mismatch');
end
clear tmp

if ~isempty(setdiff(hfa{1}.label,hfa{2}.label))
    error('ERROR: channels do not match between the two analyses!');
end

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim_S;
hfa{1} = ft_selectdata(cfg_trim,hfa{1});
cfg_trim.latency = plt_vars.plt_lim_R;
hfa{2} = ft_selectdata(cfg_trim,hfa{2});

% Compute mean RT
RT_mean = mean(round(srate*trial_info.response_time));
% Add in the baseline offset to plot correctly
RT_mean = RT_mean-plt_vars.plt_lim_S(1)*srate;

%% Plot Results
fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/actv/SR/' an_id_s '-' an_id_r '/'];
if ~exist(fig_dir,'dir')
    [~,~] = mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(hfa{1}.label)
    % Plot parameters
    fig_name = [SBJ '_actv_SR_' hfa{1}.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 0.5],'Visible',fig_vis); %twice as wide for the double plot
    plot_info.fig        = gcf;
    plot_info.x_step     = plt_vars.x_step_sz*srate;
    plot_info.legend_loc = plt_vars.legend_loc;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.sig_color  = plt_vars.sig_color;
    % Condition plotting params
    data_info.name       = {'HFA'};
    data_info.style      = {'-'};
    data_info.color      = {'k'};
    data_info.alpha      = plt_vars.errbar_alpha;
    
    for ep_ix = 1:numel(evnt_lab)
        subplot(1,numel(evnt_lab),ep_ix);
        plot_info.ax     = gca;
        plot_info.title  = [hfa{ep_ix}.label{ch_ix} ':' evnt_lab{ep_ix}];
        plot_info.legend = plt_vars.legend;
        if strcmp(evnt_lab{ep_ix},'stim')
            plot_info.x_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            % Stimulus plotting params
            event_info.name  = {evnt_lab{ep_ix}, 'Mean RT'};
            event_info.color = {[0 0 0], [0 0 0]};
            event_info.width = [plt_vars.evnt_width,plt_vars.evnt_width];
            event_info.style = {'-', plt_vars.evnt_style};
            event_info.time  = [-plt_vars.plt_lim_S(1)*srate, RT_mean];
        else
            plot_info.x_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Stimulus plotting params
            event_info.name  = {evnt_lab{ep_ix}};
            event_info.color = {plt_vars.evnt_color};
            event_info.width = plt_vars.evnt_width;
            event_info.style = {plt_vars.evnt_style};
            event_info.time  = -plt_vars.plt_lim_R(1)*srate;
        end
        
        % Compute means and variance
        avg = squeeze(mean(hfa{ep_ix}.powspctrm(:,ch_ix,:,:),1))';
        var = squeeze(std(hfa{ep_ix}.powspctrm(:,ch_ix,:,:),[],1)./sqrt(size(hfa{ep_ix}.powspctrm,1)))';
        % Find significant time periods
        if ~isempty(strmatch(hfa{ep_ix}.label{ch_ix},actv_ch{ep_ix},'exact'))
            % Find significant epoch indices
            actv_ix = strmatch(hfa{ep_ix}.label{ch_ix},actv_ch{ep_ix},'exact');
            sig_chunks = NaN(size(actv_ch_epochs{ep_ix}{actv_ix}));
            for win_ix = 1:size(actv_ch_epochs{ep_ix}{actv_ix},1)
                sig_chunks(win_ix,1) = find(hfa{ep_ix}.time==actv_ch_epochs{ep_ix}{actv_ix}(win_ix,1));
                sig_chunks(win_ix,2) = find(hfa{ep_ix}.time==actv_ch_epochs{ep_ix}{actv_ix}(win_ix,2));
            end
            
            fprintf('%s -- %i ACTIVE CLUSTERS FOUND, plotting with significance shading...\n',...
                hfa{ep_ix}.label{ch_ix},size(sig_chunks,1));
            fn_plot_ts_error_bar_sig(plot_info,avg,var,sig_chunks,event_info,data_info);
        else
            fprintf('%s -- NO ACTIVE CLUSTERS FOUND, plotting without significance shading...\n',...
                        hfa{ep_ix}.label{ch_ix});
            fn_plot_ts_error_bar(plot_info,avg,var,event_info,data_info);
        end
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

end
