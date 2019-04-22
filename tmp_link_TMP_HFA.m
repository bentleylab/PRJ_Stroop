function SBJ10b_HFA_link_sig_ANOVA(SBJs,pipeline_id,stat_id,an_id_s,an_id_r)
% creates symlinks for all of the significant files
% Load HFA analysis results for active and condition-differentiating
%   epochs, plot a summary of those time period per electrode
% clear all; %close all;
% fig_filetype = 'png';
% if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Prep variables
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id_s '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

%!!! HACK:
actv_win = 100;

% Get condition info
[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
cond_lab = {'corr(RT)', grp_lab{:}};


%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Active channels
    actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id_s,'_mn',num2str(actv_win),'.mat');
    s_actv = load(actv_filename,'actv_ch');
    actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id_r,'_mn',num2str(actv_win),'.mat');
    r_actv = load(actv_filename,'actv_ch');
    actv_ch = union(s_actv.actv_ch,r_actv.actv_ch);
    
%     actv_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/actv/' an_id_s '-' an_id_r '/sig_ch/'];
%     if ~exist(actv_dir,'dir')
%         mkdir(actv_dir);
%     end
%     for ch_ix = 1:numel(actv_ch)
%         cd(actv_dir);
%         link_cmd = ['ln -s ../' SBJ '_actv_SR_' actv_ch{ch_ix} '.png .'];
%         system(link_cmd);
%     end
    
    for ix = 1:2
        % Load data
        if ix==1
            an_id = an_id_s;
        else
            an_id = an_id_r;
        end
        load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id,'.mat'));
        
        % FDR correct pvalues for ANOVA
        qvals = NaN(size(w2.pval));
        for ch_ix = 1:numel(stat.label)
            [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
        end
        
        %% Aggregate results per ROI
        stack_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/stack_CNI/' an_id_s '-' an_id_r '/sig_ch/'];
        if ~exist(stack_dir,'dir')
            mkdir(stack_dir);
        end
        ts_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/' stat_id '/SR/' an_id_s '-' an_id_r '/sig_ch/'];
        if ~exist(ts_dir,'dir')
            mkdir(ts_dir);
        end
        for ch_ix = 1:numel(stat.label)
            % Get RT correlation onset
%             if sum(squeeze(stat.mask(ch_ix,1,:)))>0
%                 cd(stack_dir);
%                 link_cmd = ['ln -s ../' SBJ '_CNI_SR_stack_' stat.label{ch_ix} '.png .'];
%                 system(link_cmd);
%             end
            
            % Get ANOVA group onsets
            for grp_ix = 1:numel(grp_lab)
                if strcmp(grp_lab{grp_ix},'CNI') && any(squeeze(qvals(grp_ix,ch_ix,:))<0.05)
%                     sig_onsets = stat.time(win_lim(squeeze(qvals(grp_ix,ch_ix,:))<0.05,1));
%                     if strcmp(event_lab,'resp') && (sig_onsets(1)<0)
%                         all_onsets{sbj_ix,roi_ix,grp_ix+1} = [all_onsets{sbj_ix,roi_ix,grp_ix+1} sig_onsets(1)];
%                         all_onset_elec_lab{sbj_ix,roi_ix,grp_ix+1} = [all_onset_elec_lab{sbj_ix,roi_ix,grp_ix+1} stat.label(ch_ix)];
%                     elseif strcmp(event_lab,'stim') && (sig_onsets(1)<mean_RTs(sbj_ix))
%                         all_onsets{sbj_ix,roi_ix,grp_ix+1} = [all_onsets{sbj_ix,roi_ix,grp_ix+1} sig_onsets(1)];
%                         all_onset_elec_lab{sbj_ix,roi_ix,grp_ix+1} = [all_onset_elec_lab{sbj_ix,roi_ix,grp_ix+1} stat.label(ch_ix)];
%                     end
                    cd(stack_dir);
                    link_cmd = ['ln -s ../' SBJ '_CNI_SR_stack_' stat.label{ch_ix} '.png .'];
                    system(link_cmd);
%                     cd(ts_dir);
%                     link_cmd = ['ln -s ../' SBJ '_ANOVA_' stat_id '_SR_' stat.label{ch_ix} '.svg .'];
%                     system(link_cmd);
                end
            end
        end
    end
    clear SBJ SBJ_vars hfa stat w2
end
    
end
