function SBJ08b_HFA_plot_SR_stats(SBJ,conditions,an_id_s,an_id_r,plt_id,save_fig,fig_vis,fig_ftype)
% Plots HFA means with significance between conditions highlighted

%% Data Preparation
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);
event_lab = {'stim', 'resp'};

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

hfa_fname1 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_s,'.mat');
hfa_fname2 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_r,'.mat');
tmp = load(hfa_fname1,'hfa'); hfa{1} = tmp.hfa;
tmp = load(hfa_fname2,'hfa'); hfa{2} = tmp.hfa;
stats_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_s,'_',conditions,'.mat');
stats_filename2 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_r,'_',conditions,'.mat');
tmp = load(stats_filename1,'stat'); stat{1} = tmp.stat;
tmp = load(stats_filename2,'stat'); stat{2} = tmp.stat;
clear tmp

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(hfa{1}.time)-1)/(hfa{1}.time(end)-hfa{1}.time(1));
if ~isempty(setdiff(stat{1}.label,stat{2}.label))
    error('ERROR: channels do not match between the two analyses!');
end

%% Separate Trials by Condition
hfa_all = hfa;
hfa     = cell([2 numel(cond_lab)]);
for sr_ix = 1:2
    for cond_ix = 1:numel(cond_lab)
        cfgs = [];
        cfgs.trials = find(fn_condition_index(cond_lab{cond_ix}, trial_info.condition_n));
        hfa{sr_ix,cond_ix} = ft_selectdata(cfgs,hfa_all{sr_ix});
    end
end

%% Load elec structs
% % Dx for anatomy
% elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_pat_Dx.mat'];
% tmp = load(elec_atlas_fname); elec_dx = tmp.elec;
% cfg = []; cfg.channel = stat{1}.label;
% elec_dx = fn_select_elec(cfg,elec_dx);
% elec_dx = fn_reorder_elec(elec_dx,stat{1}.label);
% elec_dx.roi = fn_atlas2roi_labels(elec_dx.atlas_label,'Dx','gROI');
% % Yeo7 for network
% elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_mni_v_Yeo7.mat'];
% tmp = load(elec_atlas_fname); elec_yeo7 = tmp.elec;
% elec_yeo7 = fn_select_elec(cfg,elec_yeo7);
% elec_yeo7.roi = fn_atlas2roi_labels(elec_yeo7.atlas_label,'Yeo7','Yeo7');

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim_S;
hfa{1,1} = ft_selectdata(cfg_trim,hfa{1,1});
hfa{1,2} = ft_selectdata(cfg_trim,hfa{1,2});
cfg_trim.latency = plt_vars.plt_lim_R;
hfa{2,1} = ft_selectdata(cfg_trim,hfa{2,1});
hfa{2,2} = ft_selectdata(cfg_trim,hfa{2,2});

% Compute mean RT per condition
RTs = round(sample_rate*trial_info.response_time);
for cond_ix = 1:numel(cond_lab)
    RT_means{cond_ix} = mean(RTs(fn_condition_index(cond_lab{cond_ix}, trial_info.condition_n)==1));
    % Add in the baseline offset to plot correctly
    RT_means{cond_ix} = RT_means{cond_ix}-plt_vars.plt_lim_S(1)*sample_rate;
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/' conditions '/SR/' an_id_s '-' an_id_r '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
sig_ch = {};
for ch_ix = 1:numel(stat{1}.label)
    % Plot parameters
    fig_name = [SBJ '_' conditions '_SR_' stat{1}.label{ch_ix}];
%     [plot_rc,~] = fn_num_subplots(numel(stat.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 0.5],'Visible',fig_vis); %twice as wide for the double plot
    plot_info.fig        = gcf;
    plot_info.x_step     = plt_vars.x_step_sz*sample_rate;
    plot_info.legend_loc = plt_vars.legend_loc;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.sig_color  = plt_vars.sig_color;
    % Condition plotting params
    cond_info.name       = cond_lab;
    cond_info.style      = cond_style;
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_lab)]);
    
    for sr_ix = 1:2
        subplot(1,2,sr_ix);
        plot_info.ax     = gca;
        plot_info.title  = stat{sr_ix}.label{ch_ix};% '-(Dx,' elec_dx.roi{ch_ix} ',' ...
                                    %elec_dx.atlas_label{ch_ix} ')-(Yeo7,' elec_yeo7.roi{ch_ix} ')'];
        plot_info.legend = plt_vars.legend;
        if strcmp(event_lab{sr_ix},'stim')
            plot_info.x_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{sr_ix}, cond_lab{:}};
            event_info.color = {[0 0 0], cond_colors{:}};
            event_info.width = repmat(plt_vars.evnt_width,[1 numel(event_info.name)]);
            event_info.style = repmat({plt_vars.evnt_style},[1 numel(event_info.name)]);
            event_info.time  = [-plt_vars.plt_lim_S(1)*sample_rate, RT_means{:}];
        else
            plot_info.x_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{sr_ix}};
            event_info.width = plt_vars.evnt_width;
            event_info.color = {plt_vars.evnt_color};
            event_info.style = {plt_vars.evnt_style};
            event_info.time  = -plt_vars.plt_lim_R(1)*sample_rate;
        end
        
        % Compute means and variance
        means = NaN([numel(cond_lab) size(hfa{sr_ix,1}.powspctrm,4)]);
        var = NaN([numel(cond_lab) size(hfa{sr_ix,1}.powspctrm,4)]);
        for cond_ix = 1:numel(cond_lab)
            means(cond_ix,:) = squeeze(mean(hfa{sr_ix,cond_ix}.powspctrm(:,ch_ix,:,:),1));
            var(cond_ix,:) = squeeze(std(hfa{sr_ix,cond_ix}.powspctrm(:,ch_ix,:,:),[],1)./...
                                                    sqrt(size(hfa{sr_ix,cond_ix}.powspctrm,1)))';
        end
        % Find significant time periods
        if sum(stat{sr_ix}.mask(ch_ix,:))>0
            sig_ch = {sig_ch{:} stat{sr_ix}.label{ch_ix}};
            mask_chunks = fn_find_chunks(stat{sr_ix}.mask(ch_ix,:));
            sig_chunks = mask_chunks;
            sig_chunks(stat{sr_ix}.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
            % If stat and hfa aren't on same time axis, adjust sig_chunk indices
            if (size(stat{sr_ix}.time,2)~=size(hfa{sr_ix,1}.time,2)) || ...
                                (sum(stat{sr_ix}.time==hfa{sr_ix,1}.time)~=numel(stat{sr_ix}.time))
                for chunk_ix = 1:size(sig_chunks,1)
                    sig_chunks(chunk_ix,1) = find(hfa{sr_ix,1}.time==stat{sr_ix}.time(sig_chunks(chunk_ix,1)));
                    sig_chunks(chunk_ix,2) = find(hfa{sr_ix,1}.time==stat{sr_ix}.time(sig_chunks(chunk_ix,2)));
                end
            end
            fprintf('%s -- %i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
                stat{sr_ix}.label{ch_ix},size(sig_chunks,1));
            fn_plot_ts_error_bar_sig(plot_info,means,var,sig_chunks,event_info,cond_info);
        else
            fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat{sr_ix}.label{ch_ix});
            fn_plot_ts_error_bar(plot_info,means,var,event_info,cond_info);
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

% Save out list of channels with significant differences
sig_report_filename = [fig_dir 'sig_ch_list.txt'];
sig_report = fopen(sig_report_filename,'a');
fprintf(sig_report,'%s: %s - %s\n',datestr(datetime),an_id_s,an_id_r);
fprintf(sig_report,'%s\n',sig_ch{:});
fclose(sig_report);

end
