%% Compile Results
% Create a table of all electrodes, their ROIs, time series of significant effects
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
clear all;%close all;

SBJ        = 'IR32';
cond_name  = 'CI';              % conditions to average
analyses   = {'ERP','theta','HG'};

epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
buff_lim  = [1000, 1000];            % buffers are for cutting time series and then plotting
sig_lim   = [500 500];
% Analysis parameters
HG_type      = 'wideband';                     % only for HG: 'wideband', 'multiband'
event        = 'stim';                          % 'stim'/'resp': event to lock trials
bsln_type    = {'demean','zscore','zscore'};    % 'zscore', 'demean', 'none'
bsln_event   = 's';                    % 's'/'r': event to lock baselining
bsln_lim     = [250, -50];             % ms before and after event
smooth_it    = 1;                      % smooth the result before averaging (0/1)
smooth_freq  = 10;
trial_type   = 'datapad';               % 'nanpad', 'datapad'

% Plotting parameters
save_fig      = 1;
vis_fig       = 'on';                % 'on' or 'off' to determine if plot is shown
plot_lim      = [500 500];            % ms to plot before the event and **after RT**
x_step        = 250;                  % step of x tick marks
y_step        = 5;                  % increments of increase on y axis (will round up)
event_ln_width= 2;
fig_type      = 'png';

data_id = strcat(SBJ,'_IH_CA');
SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
fig_dir = strcat('/home/knight/hoycw/PRJ_Stroop/results/multifreq/',SBJ,'/',cond_name,'/',event,'/');
all_analysis_id = '';
all_bsln_type = '';
for an_ix = 1:length(analyses)
    [analysis_id{an_ix}, env_it{an_ix}, filt_lim{an_ix}] =...
        fn_B02_analysis_params(analyses{an_ix},HG_type);
    all_analysis_id = strcat(all_analysis_id, analysis_id{an_ix}(1));
    all_bsln_type = strcat(all_bsln_type, bsln_type{an_ix}(1));
end
bsln_id  = ['_' all_bsln_type '.' bsln_event num2str(bsln_lim(1)) '.' num2str(bsln_lim(2))];
epoch_id = [num2str(floor(epoch_lim(1)),'%+d') '.' num2str(floor(epoch_lim(2)),'%+d')];
sig_id = ['sig' num2str(sig_lim(1)) '.' num2str(sig_lim(2))];
if strcmp(trial_type,'nanpad')
    trl_id = 'nan';
elseif strcmp(trial_type,'datapad')
    trl_id = 'data';
else
    error(strcat('Unknown trial averaging type: ',trial_type));
end
if smooth_it == 0
    smooth_id = '';
elseif smooth_it == 1
    smooth_id = strcat('_sm',num2str(smooth_freq));
else
    error('smooth_it not in [0,1]');
end
plot_lim_id = strcat('_xlim',num2str(plot_lim(1)),'.',num2str(plot_lim(2)));

%% Load data
% Significance times
sig_times_file = [SBJ_dir '04_proc/' data_id '_' all_analysis_id...
    '_' event '_' cond_name '_Bob.ep' epoch_id '_' trl_id '_bsln' bsln_id smooth_id '_' sig_id '_times.mat'];
%'sig_times','sig_wins','sig_pvals','ch_ids','header_ecog'
% sig_times{analysis_ix}{elec_ix} = [array of time points as double]
load(sig_times_file); 

%% Load RT info
load(strcat(SBJ_dir,'03_events/tmp_',SBJ,'_B02_max_RT.mat'));

%% Load ROI and GM/WM info
roi_filename = strcat(SBJ_dir,'04_proc/',SBJ,'_sfn_elec_ROIs.csv');
einfo = textread(roi_filename,'%s','delimiter',',');
einfo = reshape(einfo,[3, length(header_ecog.channel_labels)]); % elec_name,ROI,tissue(GG,GW,GW,WW,BS)
gm_elec_idx = logical(zeros([1 size(einfo,2)]));
for e_ix = 1:size(einfo,2)
    if ~isempty(strfind(einfo{3,e_ix},'G'))
        gm_elec_idx(e_ix) = true;
    end
end
einfo_gm = einfo(:,gm_elec_idx);
roi_list = unique(einfo_gm(2,:));
% roi_elecs = {};
% for e_ix = 1:size(einfo,2)
%     roi_elecs{roi_ix} = 
% end
%% Plot significance time points for each region and analysis
fig_file = [data_id '_allE_sig_ts_' all_analysis_id...
    '_' event '_' cond_name '_Bob.ep' epoch_id '_' trl_id...
    '_bsln' bsln_id smooth_id plot_lim_id];
fig_height = 1;%length(analyses)/3;
figure('Name',fig_file,'units','normalized',...
    'outerposition',[0 0 1 fig_height],'Visible',vis_fig);

if strcmp(event,'stim')
    trial_len = buff_lim(1)+max_RT+buff_lim(2)+1;
    win_off = buff_lim(1)+max_RT+plot_lim(2)+1;         % Capture even longest trials
else
    trial_len = buff_lim(1)+buff_lim(2)+1;
    win_off = buff_lim(1)+plot_lim(2)+1;
end
win_on  = buff_lim(1)-plot_lim(1)+1;
x_lab   = -plot_lim(1):x_step:max_RT+plot_lim(2);
y_max   = 0;
for an_ix = 1:length(analyses)
    for roi_ix = 1:length(roi_list)
        roi_idx = ~cellfun('isempty',strfind(einfo(2,:),roi_list{roi_ix}));
        roi_idx(~gm_elec_idx) = 0;
        if sum([sig_times{an_ix}{roi_idx}])>0
            [overlap_cnt,~] = histcounts([sig_times{an_ix}{roi_idx}],unique([sig_times{an_ix}{roi_idx}]));
            y_max = max([y_max max(overlap_cnt)]);
        end
    end
end
y_max = y_step*floor(y_max/y_step)+y_step;

for an_ix = 1:length(analyses)
    for roi_ix = 1:length(roi_list)
        subplot(length(analyses),length(roi_list),...
            fn_RC2subplot_ix(length(analyses),length(roi_list),an_ix,roi_ix));
        sig_ts = zeros([1 trial_len]);
        n_names = 0;
        elec_names = '';
        for e_ix = 1:size(einfo,2)
            % must be (1) in the ROI (2) have GM and (3) be significant
            if (strcmp(einfo{2,e_ix},roi_list{roi_ix})) && (gm_elec_idx(e_ix)==1) && ~isempty((sig_times{an_ix}{e_ix}))
                sig_ts(sig_times{an_ix}{e_ix}) = sig_ts(sig_times{an_ix}{e_ix})+1;
                n_names = n_names+1;
                % Shorten the name
                [~, uniq_pos, ~] = unique(einfo{1,e_ix});
                name = einfo{1,e_ix}(sort(uniq_pos));
                if n_names==1
                    elec_names = name;
                else
                    elec_names = [elec_names ', ' name];
%                     if mod(n_names,4)==0
%                         elec_names = [elec_names '\n'];
%                     end
                end
            end
        end
        plot(sig_ts(win_on:win_off));
        ax = gca;
        ax.YLim = [0 y_max];
        event_line = line([plot_lim(1) plot_lim(1)],ylim,...
            'LineWidth',event_ln_width,'Color','k');
        % Axes and Labels
        ax.XLim = [0,win_off-win_on];
        ax.XTick = 0:x_step:win_off;
        ax.XTickLabel = x_lab;
%         main_lines = [main_lines event_line];
%         if an_ix==1
%             legend(main_lines,cond_lab_legend{:},event,'Location','southeast');%,'inc-con'
%         end
        title(sprintf('%s - %s:\n%s',analyses{an_ix},roi_list{roi_ix},elec_names));
    end
end

%% Save figure
fig_filename = strcat(fig_dir,fig_file);
if save_fig ==1
    fprintf('Saving %s\n',fig_filename);
    eval(['export_fig ' fig_filename]);
end



% %% Write to csv
% fid = fopen(roi_filename, 'w');
% for lab_ix = 1:length(header_ecog.channel_labels)
%     fprintf(fid, '%s\n', header_ecog.channel_labels{lab_ix});
% end
% fclose(fid);
