function SBJ08b_HFA_plot_SR_stats_svg(SBJ,elec,conditions,pipeline_id,an_id_s,an_id_r,plt_id,save_fig,fig_vis)
% Plots ERPs computed in SBJ07a_ERP_stats
% clear all; %close all;

fig_filetype = 'svg';
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

einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' pipeline_id '.mat'];
load(einfo_filename);
% Electrode Info Table:
%   label- name of electrode
%   ROI- specific region
%   gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
%   ROI2- specific region of second electrode
%   tissue- primary tissue type
%   GM weight- percentage of electrode pair in GM
%   Out- 0/1 flag for whether this is partially out of the brain
einfo_ix = strmatch(elec,einfo(:,1),'exact');
roi = einfo{einfo_ix,2};

stats_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id_s,'.mat');
stats_filename2 = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id_r,'.mat');
tmp = load(stats_filename1,'stat'); stat{1} = tmp.stat;
tmp = load(stats_filename2,'stat'); stat{2} = tmp.stat;
tmp = load(stats_filename1,'hfa'); hfa{1,1} = tmp.hfa{1}; hfa{1,2} = tmp.hfa{2};
tmp = load(stats_filename2,'hfa'); hfa{2,1} = tmp.hfa{1}; hfa{2,2} = tmp.hfa{2};
clear tmp

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(hfa{1,1}.time)-1)/(hfa{1,1}.time(end)-hfa{1,1}.time(1));
if ~isempty(setdiff(stat{1}.label,stat{2}.label))
    error('ERROR: channels do not match between the two analyses!');
end

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.channel = elec;
stat{1} = ft_selectdata(cfg_trim,stat{1});
stat{2} = ft_selectdata(cfg_trim,stat{2});
cfg_trim.latency = plt_vars.plt_lim_S;
hfa{1,1} = ft_selectdata(cfg_trim,hfa{1,1});
hfa{1,2} = ft_selectdata(cfg_trim,hfa{1,2});
cfg_trim.latency = plt_vars.plt_lim_R;
hfa{2,1} = ft_selectdata(cfg_trim,hfa{2,1});
hfa{2,2} = ft_selectdata(cfg_trim,hfa{2,2});

% Compute mean RT per condition
RTs = round(1000*trial_info.response_time); % converts sec to ms
for cond_ix = 1:numel(cond_lab)
    RT_means{cond_ix} = mean(RTs(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1));
    % Add in the baseline offset to plot correctly
    RT_means{cond_ix} = RT_means{cond_ix}-plt_vars.plt_lim_S(1)*1000;
end

%% Plot Results
fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/' conditions '/SR/' an_id_s '-' an_id_r '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
% Plot parameters
fig_name = strcat(SBJ,'_',conditions,'_SR_',roi,'_',elec);
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.4 0.5],'Visible',fig_vis); %twice as wide for the double plot
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

for set_ix = 1:2
    subplot(1,2,set_ix);
    plot_info.ax     = gca;
    plot_info.title  = [event_lab{set_ix} '-Locked'];
    plot_info.title(1) = upper(plot_info.title(1));
        if set_ix==2
            plot_info.legend = 1;
        else
            plot_info.legend = 2;
        end
        if strcmp(event_lab{set_ix},'stim')
            plot_info.x_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{set_ix}, cond_lab{:}};
            event_info.color = {[0 0 0], cond_colors{:}};
            event_info.width = repmat(plt_vars.evnt_width,[1 numel(event_info.name)]);
            event_info.style = repmat({plt_vars.evnt_style},[1 numel(event_info.name)]);
            event_info.style{1} = '-';
            event_info.time  = [-plt_vars.plt_lim_S(1)*sample_rate, RT_means{:}];
        else
            plot_info.x_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{set_ix}};
            event_info.width = plt_vars.evnt_width;
            event_info.color = {plt_vars.evnt_color};
            event_info.style = {plt_vars.evnt_style};
            event_info.time  = -plt_vars.plt_lim_R(1)*sample_rate;
        end
        
        % Compute means and variance
        means = NaN([numel(cond_lab) size(hfa{set_ix,1}.powspctrm,4)]);
        var = NaN([numel(cond_lab) size(hfa{set_ix,1}.powspctrm,4)]);
        for cond_ix = 1:numel(cond_lab)
            means(cond_ix,:) = squeeze(mean(hfa{set_ix,cond_ix}.powspctrm(:,1,:,:),1));
            var(cond_ix,:) = squeeze(std(hfa{set_ix,cond_ix}.powspctrm(:,1,:,:),[],1)./...
                sqrt(size(hfa{set_ix,cond_ix}.powspctrm,1)))';
        end
        % Find significant time periods
        if sum(stat{set_ix}.mask(1,:))>0
            mask_chunks = fn_find_chunks(stat{set_ix}.mask(1,:));
            sig_chunks = mask_chunks;
            sig_chunks(stat{set_ix}.mask(1,sig_chunks(:,1))==0,:) = [];
            % If stat and hfa aren't on same time axis, adjust sig_chunk indices
            if (size(stat{set_ix}.time,2)~=size(hfa{set_ix,1}.time,2)) || ...
                    (sum(stat{set_ix}.time==hfa{set_ix,1}.time)~=numel(stat{set_ix}.time))
                for chunk_ix = 1:size(sig_chunks,1)
                    sig_chunks(chunk_ix,1) = find(hfa{set_ix,1}.time==stat{set_ix}.time(sig_chunks(chunk_ix,1)));
                    sig_chunks(chunk_ix,2) = find(hfa{set_ix,1}.time==stat{set_ix}.time(sig_chunks(chunk_ix,2)));
                end
            end
            fprintf('%s -- %i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
                elec,size(sig_chunks,1));
            fn_plot_ts_error_bar_sig(plot_info,means,var,sig_chunks,event_info,cond_info);
        else
            fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',elec);
            fn_plot_ts_error_bar(plot_info,means,var,event_info,cond_info);
        end
end

%% Save figure
if save_fig
    fig_filename = [fig_dir fig_name '.' fig_filetype];
    fprintf('Saving %s\n',fig_filename,'svg');
    saveas(gcf,fig_filename);
    %eval(['export_fig ' fig_filename]);
end


end
