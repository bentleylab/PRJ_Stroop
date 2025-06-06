function SBJ07b_ERP_plot_SR_but_stats(SBJ,conditions,an_id_s,an_id_r,plt_id,save_fig,fig_vis)
% Plots ERPs computed in SBJ07a_ERP_stats
%   Butterfly plots per condition with mean overlay, and normal condition
%   means overlay in a separate se tof subplots
% clear all; %close all;

fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

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

[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);
event_lab = {'stim', 'resp'};
% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

stats_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id_s,'.mat');
stats_filename2 = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id_r,'.mat');
tmp = load(stats_filename1,'stat'); stat{1} = tmp.stat;
tmp = load(stats_filename2,'stat'); stat{2} = tmp.stat;
tmp = load(stats_filename1,'roi_erp'); erp{1,1} = tmp.roi_erp{1}; erp{1,2} = tmp.roi_erp{2};
tmp = load(stats_filename2,'roi_erp'); erp{2,1} = tmp.roi_erp{1}; erp{2,2} = tmp.roi_erp{2};
clear tmp

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(erp{1,1}.time)-1)/(erp{1,1}.time(end)-erp{1,1}.time(1));
if ~isempty(setdiff(stat{1}.label,stat{2}.label))
    error('ERROR: channels do not match between the two analyses!');
end

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim_S;
erp{1,1} = ft_selectdata(cfg_trim,erp{1,1});
erp{1,2} = ft_selectdata(cfg_trim,erp{1,2});
cfg_trim.latency = plt_vars.plt_lim_R;
erp{2,1} = ft_selectdata(cfg_trim,erp{2,1});
erp{2,2} = ft_selectdata(cfg_trim,erp{2,2});

% Compute mean RT per condition
RTs = round(1000*trial_info.response_time); % converts sec to ms
for cond_ix = 1:numel(cond_lab)
    RT_means{cond_ix} = mean(RTs(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1));
    % Add in the baseline offset to plot correctly
    RT_means{cond_ix} = RT_means{cond_ix}-plt_vars.plt_lim_S(1)*1000;
end

%% Plot Results
% NO! I will plot whatever channels I ran the stats on (why else did I run them?)
% % Select data to plot only ROI channels
% cfgs = [];
% cfgs.channel = SBJ_vars.ch_lab.ROI;
% stat = ft_selectdata(cfgs,stat);
% for an_ix = 1:numel(cond_lab)
%     erp{an_ix} = ft_selectdata(cfgs,erp{an_ix});
% end

fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/ERP/' SBJ '/' conditions '/SR_trl/' an_id_s '-' an_id_r '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
sig_ch = {};
for ch_ix = 1:numel(stat{1}.label)
    % Plot parameters
%     probe_name = stat.label{ch_ix}(regexp(stat.label{ch_ix},'\D'));
%     probe_name(strfind(probe_name,'-')) = [];
    
    fig_name = [SBJ '_' conditions '_SR_but_' stat{1}.label{ch_ix}];
%     [plot_rc,~] = fn_num_subplots(numel(stat.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);   %this size is for single plots
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
    cond_info.name       = cond_lab;
    cond_info.style      = cond_style;
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_lab)]);
    
    % Butterfly trial plots
    for epoch_ix = 1:2
        for cond_ix = 1:numel(cond_lab)
            subplot(numel(cond_lab)+1,2,fn_rowcol2subplot_ix(numel(cond_lab)+1,2,cond_ix,epoch_ix));
            plot_info.ax     = gca;
            plot_info.title  = [stat{epoch_ix}.label{ch_ix} ': ' cond_lab{cond_ix} '-' event_lab{epoch_ix}];
            plot_info.legend = plt_vars.legend;
            if strcmp(event_lab{epoch_ix},'stim')
                plot_info.x_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
                % Stimulus plotting params
                event_info.name  = {event_lab{epoch_ix}, cond_lab{cond_ix}};
                event_info.color = {[0 0 0], cond_colors{cond_ix}};
                event_info.width = repmat(plt_vars.evnt_width,[1 numel(event_info.name)]);
                event_info.style = repmat({plt_vars.evnt_style},[1 numel(event_info.name)]);
                event_info.time  = [-plt_vars.plt_lim_S(1)*sample_rate, RT_means{cond_ix}];
            else
                plot_info.x_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
                % Stimulus plotting params
                event_info.name  = {event_lab{epoch_ix}};
                event_info.width = plt_vars.evnt_width;
                event_info.color = {plt_vars.evnt_color};
                event_info.style = {plt_vars.evnt_style};
                event_info.time  = -plt_vars.plt_lim_R(1)*sample_rate;
            end
            
            % Grab trials, means, and variance
            trials = squeeze(erp{epoch_ix,cond_ix}.trial(:,ch_ix,:));
            fn_plot_ts_but(plot_info,trials,event_info,cond_info,cond_ix);
            clear event_info
        end
    end
    
    % Mean error bar plots
    for epoch_ix = 1:2
        subplot(numel(cond_lab)+1,2,fn_rowcol2subplot_ix(numel(cond_lab)+1,2,numel(cond_lab)+1,epoch_ix));
        plot_info.ax     = gca;
        plot_info.title  = [stat{epoch_ix}.label{ch_ix} ':' event_lab{epoch_ix}];
        plot_info.legend = plt_vars.legend;
        if strcmp(event_lab{epoch_ix},'stim')
            plot_info.x_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{epoch_ix}, cond_lab{:}};
            event_info.color = {[0 0 0], cond_colors{:}};
            event_info.width = repmat(plt_vars.evnt_width,[1 numel(event_info.name)]);
            event_info.style = repmat({plt_vars.evnt_style},[1 numel(event_info.name)]);
            event_info.time  = [-plt_vars.plt_lim_S(1)*sample_rate, RT_means{:}];
        else
            plot_info.x_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{epoch_ix}};
            event_info.width = plt_vars.evnt_width;
            event_info.color = {plt_vars.evnt_color};
            event_info.style = {plt_vars.evnt_style};
            event_info.time  = -plt_vars.plt_lim_R(1)*sample_rate;
        end
        
        % Grab means and variance
        means  = NaN([numel(cond_lab) size(erp{epoch_ix,1}.avg,2)]);
        var    = NaN([numel(cond_lab) size(erp{epoch_ix,1}.avg,2)]);
        for an_ix = 1:numel(cond_lab)
            means(an_ix,:) = erp{epoch_ix,an_ix}.avg(ch_ix,:);
            var(an_ix,:) = squeeze(std(erp{epoch_ix,an_ix}.trial(:,ch_ix,:),[],1)./...
                                                sqrt(size(erp{epoch_ix,an_ix}.trial,1)))';
        end
        
        % Find significant time periods
        if sum(stat{epoch_ix}.mask(ch_ix,:))>0
            sig_ch = {sig_ch{:} stat{epoch_ix}.label{ch_ix}};
            mask_chunks = fn_find_chunks(stat{epoch_ix}.mask(ch_ix,:));
            sig_chunks = mask_chunks;
            sig_chunks(stat{epoch_ix}.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
            % If stat and erp aren't on same time axis, adjust sig_chunk indices
            if (size(stat{epoch_ix}.time,2)~=size(erp{epoch_ix,1}.time,2)) || ...
                                (sum(stat{epoch_ix}.time==erp{epoch_ix,1}.time)~=numel(stat{epoch_ix}.time))
                for chunk_ix = 1:size(sig_chunks,1)
                    sig_chunks(chunk_ix,1) = find(erp{epoch_ix,1}.time==stat{epoch_ix}.time(sig_chunks(chunk_ix,1)));
                    sig_chunks(chunk_ix,2) = find(erp{epoch_ix,1}.time==stat{epoch_ix}.time(sig_chunks(chunk_ix,2)));
                end
            end
            fprintf('%s -- %i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
                stat{epoch_ix}.label{ch_ix},size(sig_chunks,1));
            fn_plot_ts_error_bar_sig(plot_info,means,var,sig_chunks,event_info,cond_info);
        else
            fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat{epoch_ix}.label{ch_ix});
            fn_plot_ts_error_bar(plot_info,means,var,event_info,cond_info);
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

% Save out list of channels with significant differences
sig_report_filename = [fig_dir 'sig_ch_list.txt'];
sig_report = fopen(sig_report_filename,'a');
fprintf(sig_report,'%s - %s\n',an_id_s,an_id_r);
fprintf(sig_report,'%s\n',sig_ch{:});
fclose(sig_report);

end