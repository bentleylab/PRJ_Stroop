function SBJ08b_HFA_plot_SR_but_actv(SBJ,an_id_s,an_id_r,actv_win,plt_id,save_fig,fig_vis)
% Plots both stimulus- and response-locked HFA computed in SBJ08a_HFA_actv
% clear all; %close all;

fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

event_lab = {'stim', 'resp'};
% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

hfa_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id_s,'_mn',actv_win,'.mat');
hfa_filename2 = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id_r,'_mn',actv_win,'.mat');
tmp = load(hfa_filename1,'hfa'); hfa{1} = tmp.hfa;
tmp = load(hfa_filename2,'hfa'); hfa{2} = tmp.hfa;
tmp = load(hfa_filename1,'actv_ch'); actv_ch{1} = tmp.actv_ch;
tmp = load(hfa_filename2,'actv_ch'); actv_ch{2} = tmp.actv_ch;
tmp = load(hfa_filename1,'actv_ch_epochs'); actv_ch_epochs{1} = tmp.actv_ch_epochs;
tmp = load(hfa_filename2,'actv_ch_epochs'); actv_ch_epochs{2} = tmp.actv_ch_epochs;
clear tmp

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(hfa{1}.time)-1)/(hfa{1}.time(end)-hfa{1}.time(1));
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
RT_mean = mean(round(1000*trial_info.response_time)); % converts sec to ms
% Add in the baseline offset to plot correctly
RT_mean = RT_mean-plt_vars.plt_lim_S(1)*1000;

%% Plot Results
fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/actv/SR_but/' an_id_s '-' an_id_r '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(hfa{1}.label)
    % Plot parameters
    fig_name = [SBJ '_actv_SR_but_' hfa{1}.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    plot_info.fig        = gcf;
    plot_info.x_step     = plt_vars.x_step_sz*sample_rate;
    plot_info.legend_loc = plt_vars.legend_loc;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.sig_color  = plt_vars.sig_color;
    plot_info.but_width  = plt_vars.but_width;
    plot_info.but_color  = plt_vars.but_color;
    plot_info.but_alpha  = plt_vars.but_alpha;
    plot_info.mean_width = plt_vars.mean_width;
    % Condition plotting params
    data_info.name       = {'HFA'};
    data_info.style      = {'-'};
    data_info.color      = {'b'};
    data_info.alpha      = plt_vars.errbar_alpha;
    
    for ep_ix = 1:numel(event_lab)
        subplot(2,numel(event_lab),fn_rowcol2subplot_ix(2,2,1,ep_ix));
        plot_info.ax     = gca;
        plot_info.title  = [hfa{ep_ix}.label{ch_ix} ':' event_lab{ep_ix} ' trials'];
        plot_info.legend = plt_vars.legend;
        if strcmp(event_lab{ep_ix},'stim')
            plot_info.x_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{ep_ix}, 'Mean RT'};
            event_info.color = {[0 0 0], [0 0 0]};
            event_info.width = [plt_vars.evnt_width,plt_vars.evnt_width];
            event_info.style = {'-', plt_vars.evnt_style};
            event_info.time  = [-plt_vars.plt_lim_S(1)*sample_rate, RT_mean];
        else
            plot_info.x_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{ep_ix}};
            event_info.color = {plt_vars.evnt_color};
            event_info.width = plt_vars.evnt_width;
            event_info.style = {plt_vars.evnt_style};
            event_info.time  = -plt_vars.plt_lim_R(1)*sample_rate;
        end
        
        % Grab trials and plot them
        trials = squeeze(hfa{ep_ix}.powspctrm(:,ch_ix,:,:));
        fn_plot_ts_but(plot_info,trials,event_info,data_info,1);
        
        % Now compute and plot means by themselves
        subplot(2,numel(event_lab),fn_rowcol2subplot_ix(2,2,2,ep_ix));
        plot_info.ax     = gca;
        plot_info.title  = [hfa{ep_ix}.label{ch_ix} ':' event_lab{ep_ix} ' means'];
        
        % Compute means and variance
        avg = squeeze(mean(hfa{ep_ix}.powspctrm(:,ch_ix,:,:),1))';
        var   = squeeze(std(hfa{ep_ix}.powspctrm(:,ch_ix,:,:),[],1)./sqrt(size(hfa{ep_ix}.powspctrm,1)))';
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
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

end
