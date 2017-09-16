function SBJ08c_HFA_plot_SR_stats(SBJ,conditions,an_id1,an_id2,plt_id,save_fig,fig_vis)
% Plots ERPs computed in SBJ07a_ERP_stats
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
an_vars_cmd1 = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id1 '_vars.m'];
eval(an_vars_cmd1);
plt_lim1 = plt_lim;
evnt_type1 = event_type;
an_vars_cmd2 = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id2 '_vars.m'];
eval(an_vars_cmd2);
plt_lim2 = plt_lim;
evnt_type2 = event_type;

[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);

stats_filename1 = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_ROI_',an_id1,'.mat');
data1 = load(stats_filename1);
stats_filename2 = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_ROI_',an_id2,'.mat');
data2 = load(stats_filename2);
%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(data1.roi_hfa{1}.time)-1)/(data1.roi_hfa{1}.time(end)-data1.roi_hfa{1}.time(1));
if ~isempty(setdiff(data1.stat.label,data2.stat.label))
    error('ERROR: channels do not match between data1 and data2');
end

%% Plot Results
% NO! I will plot whatever channels I ran the stats on (why else did I run them?)
% % Select data to plot only ROI channels
% cfgs = [];
% cfgs.channel = SBJ_vars.ch_lab.ROI;
% stat = ft_selectdata(cfgs,stat);
% for an_ix = 1:numel(cond_lab)
%     roi_hfa{an_ix} = ft_selectdata(cfgs,roi_hfa{an_ix});
% end

fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
sig_ch = {};
for ch_ix = 1:numel(stat.label)
    % Plot parameters
%     probe_name = stat.label{ch_ix}(regexp(stat.label{ch_ix},'\D'));
%     probe_name(strfind(probe_name,'-')) = [];
    
    fig_name = [SBJ '_' conditions '_SR_' an_id '_' stat.label{ch_ix}];
%     [plot_rc,~] = fn_num_subplots(numel(stat.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
    plot_info.fig        = gcf;
    plot_info.x_step     = plt_vars.x_step_sz*sample_rate;
    plot_info.x_lab      = plt_lim(1):plt_vars.x_step_sz:plt_lim(2);
    plot_info.legend_loc = plt_vars.legend_loc;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.sig_color  = plt_vars.sig_color;
    % Stimulus plotting params
    event_info.time      = [-plt_lim1(1)*sample_rate -plt_lim2(1)*sample_rate];
    event_info.name      = {evnt_type1 evnt_type2};
    event_info.width     = plt_vars.evnt_width;
    event_info.color     = plt_vars.evnt_color;
    event_info.style     = plt_vars.evnt_style;
    % Condition plotting params
    cond_info.name       = cond_lab;
    cond_info.style      = cond_style;
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_lab)]);
    
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    plot_info.ax     = gca;
    plot_info.title  = stat.label{ch_ix};
    plot_info.legend = plt_vars.legend;
    
    % Compute means and variance
    means1 = NaN([numel(cond_lab) size(roi_hfa{1}.powspctrm,4)]);
    var1 = NaN([numel(cond_lab) size(roi_hfa{1}.powspctrm,4)]);
    for cond_ix = 1:numel(cond_lab)
        means1(cond_ix,:) = squeeze(mean(roi_hfa{cond_ix}.powspctrm(:,ch_ix,:,:),1));
        var1(cond_ix,:) = squeeze(std(roi_hfa{cond_ix}.powspctrm(:,ch_ix,:,:),[],1)./sqrt(size(roi_hfa{cond_ix}.powspctrm,1)))';
    end
    % Find significant time periods
    if sum(stat.mask(ch_ix,:))>0
        sig_ch = {sig_ch{:} stat.label{ch_ix}};
        mask_chunks = fn_find_chunks(stat.mask(ch_ix,:));
        sig_chunks = mask_chunks;
        sig_chunks(stat.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
        % If stat and roi_hfa aren't on same time axis, adjust sig_chunk indices
        if (size(stat.time,2)~=size(roi_hfa{1}.time,2)) || (sum(stat.time==roi_hfa{1}.time)~=numel(stat.time))
            for chunk_ix = 1:size(sig_chunks,1)
                sig_chunks(chunk_ix,1) = find(roi_hfa{1}.time==stat.time(sig_chunks(chunk_ix,1)));
                sig_chunks(chunk_ix,2) = find(roi_hfa{1}.time==stat.time(sig_chunks(chunk_ix,2)));
            end
        end
        fprintf('%s -- %i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
                                                                stat.label{ch_ix},size(sig_chunks,1));
        fn_plot_ts_error_bar_sig(plot_info,means1,var1,sig_chunks,event_info,cond_info);
    else
        fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat.label{ch_ix});
        fn_plot_ts_error_bar(plot_info,means1,var1,event_info,cond_info);
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
fprintf(sig_report,'%s\n',an_id);
fprintf(sig_report,'%s\n',sig_ch{:});
fclose(sig_report);

end