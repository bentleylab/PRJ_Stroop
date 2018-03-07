function SBJ08b_HFA_plot_actv(SBJ,an_id,actv_win,plt_id,save_fig,fig_vis)
% Plots HFA computed in SBJ08a_HFA_actv
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
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

hfa_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id,'_mn',actv_win,'.mat');
load(hfa_filename);

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(hfa.time)-1)/(hfa.time(end)-hfa.time(1));

% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim;
hfa = ft_selectdata(cfg_trim,hfa);

%% Plot Results
% NO! I will plot whatever channels I ran the stats on (why else did I run them?)
% % Select data to plot only ROI channels
% cfgs = [];
% cfgs.channel = {'LPC*'};%SBJ_vars.ch_lab.ROI;
% hfa = ft_selectdata(cfgs,hfa);

fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/actv/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(hfa.label)
    % Plot parameters
%     probe_name = hfa.label{ch_ix}(regexp(hfa.label{ch_ix},'\D'));
%     probe_name(strfind(probe_name,'-')) = [];
    
    fig_name = [SBJ '_actv_' an_id '_' hfa.label{ch_ix}];
%     [plot_rc,~] = fn_num_subplots(numel(hfa.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
    plot_info.fig        = gcf;
    plot_info.x_step     = plt_vars.x_step_sz*sample_rate;
    plot_info.x_lab      = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    plot_info.legend_loc = plt_vars.legend_loc;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.sig_color  = plt_vars.sig_color;
    % Stimulus plotting params
    event_info.time      = -plt_vars.plt_lim(1)*sample_rate;
    event_info.name      = {event_type};
    event_info.width     = plt_vars.evnt_width;
    event_info.color     = {plt_vars.evnt_color};
    event_info.style     = {plt_vars.evnt_style};
    % Condition plotting params
    data_info.name       = {'HFA'};
    data_info.style      = {'-'};
    data_info.color      = {'k'};
    data_info.alpha      = plt_vars.errbar_alpha;
    
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    plot_info.ax     = gca;
    plot_info.title  = hfa.label{ch_ix};
    plot_info.legend = plt_vars.legend;
    
    % Compute means and variance
    means = squeeze(mean(hfa.powspctrm(:,ch_ix,:,:),1))';
    var   = squeeze(std(hfa.powspctrm(:,ch_ix,:,:),[],1)./sqrt(size(hfa.powspctrm,1)))';
    
    % Find significant time periods
    if ~isempty(strmatch(hfa.label{ch_ix},actv_ch,'exact'))
        % Find significant epoch indices
        actv_ix = strmatch(hfa.label{ch_ix},actv_ch,'exact');
        sig_chunks = NaN(size(actv_ch_epochs{actv_ix}));
        for ep_ix = 1:size(actv_ch_epochs{actv_ix},1)
            sig_chunks(ep_ix,1) = find(hfa.time==actv_ch_epochs{actv_ix}(ep_ix,1));
            sig_chunks(ep_ix,2) = find(hfa.time==actv_ch_epochs{actv_ix}(ep_ix,2));
        end
        
        fprintf('%s -- %i ACTIVE CLUSTERS FOUND, plotting with significance shading...\n',...
                                                                hfa.label{ch_ix},size(sig_chunks,1));
        fn_plot_ts_error_bar_sig(plot_info,means,var,sig_chunks,event_info,data_info);
    else
        fprintf('%s -- NO ACTIVE CLUSTERS FOUND, plotting without significance shading...\n',hfa.label{ch_ix});
        fn_plot_ts_error_bar(plot_info,means,var,event_info,data_info);
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

% % Save out list of channels with significant differences
% sig_report_filename = [fig_dir 'sig_ch_list.txt'];
% sig_report = fopen(sig_report_filename,'a');
% fprintf(sig_report,'%s\n',an_id);
% fprintf(sig_report,'%s\n',sig_ch{:});
% fclose(sig_report);

end
