function SBJ08c_HFA_report_actv(SBJ,an_id,actv_win)
error('just use the SBJ08c_HFA_summary_plot_actv_cond instead, which will be a better overview');
% Load HFA analysis, write .txt report on the differences from baseline
%   (no condition differences)
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

if isnumeric(actv_win); actv_win = num2str(actv_win); end
load(strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id_s,'_mn',actv_win,'.mat'));

%% Check sign of significant epochs
actv_ch_epoch_signs = {};
for ch_ix = 1:numel(actv_ch)
    for ep_ix = 1:size(actv_ch_epochs{ch_ix},1)
        
    end
end

%% Print reports of significant channels, their signs, and their epochs
report_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/actv/'];
if ~exist(report_dir,'dir')
    mkdir(report_dir);
end

% Save out list of channels with significant activations from baseline
report_filename = [report_dir an_id '_ch_list.csv'];
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
