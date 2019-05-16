function SBJ10c_HFA_GRP_elec_cnt_tis_RT_ANOVA(SBJs,stat_id,an_id,roi_id,atlas_id,gm_thresh,pipeline_id,save_out)
% Load HFA RT correlation+ANOVA results and compare values to tissue percentages
%   plots histogram of significant electrode tissue percentages
%   prints report of # and % of sig elecs for GM and WM based on tis_thresh 
% INPUTS:
%   atlas_id [str] - name of the atlas used for ROI assignment
%       !!NOTE: automatically assumes patient space since that will be most accurate!
%   roi_id [str] - ROI groupings to convert atlas labels into
%       {'gROI','ROI'}, not any others for now...
%   tis_thresh [float] - percentage of GM necessary to be considered a valid electrode (0-1)

% clear all; %close all;
if ischar(save_out); save_out = str2num(save_out); end

%% Data Preparation
% Directories
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Prep variables
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);

% Get stat factors
[grp_lab, ~, ~] = fn_group_label_styles(model_lab);

% Get all possible ROI labels (not excluding any weird cases)
tsv_filename = [root_dir 'PRJ_Stroop/data/atlases/atlas_mappings/atlas_ROI_mappings_' atlas_id '.tsv'];
roi_file = fopen(tsv_filename, 'r');
roi_map = textscan(roi_file, '%s %s %s %s', 'HeaderLines', 1,...
    'Delimiter', '\t', 'MultipleDelimsAsOne', 0); 
fclose(roi_file);
if strcmp(roi_id,'gROI')
    map_ix = 2;
elseif strcmp(roi_id,'ROI')
    map_ix = 3;
else
    error(['roi_style unknown, must be gROI or ROI: ' roi_style]);
end
roi_list = unique(roi_map{map_ix});

% Set up electrode counts
gm_perc  = cell([numel(SBJs) 1]);
gm_bin   = cell([numel(SBJs) 1]);
grp_sig  = cell([numel(SBJs) 1]);
atlas_elec  = cell([numel(SBJs) 1]);
tis_elec  = cell([numel(SBJs) 1]);

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id,'.mat'),'stat','w2');
    
    %% FDR correct pvalues for ANOVA
    qvals = NaN(size(w2.pval));
    for ch_ix = 1:numel(stat.label)
        [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
    end
    
    %% Load ROI and GM/WM info
    elec_atlas_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_pat_' atlas_id '.mat'];
    tmp = load(elec_atlas_fname); atlas_elec{sbj_ix} = tmp.elec;
    elec_tis_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_pat_' atlas_id '_tis.mat'];
    tmp = load(elec_tis_fname); tis_elec{sbj_ix} = tmp.elec; clear tmp;
    
    % Sort elecs by stat labels
    cfgs = []; cfgs.channel = stat.label;
    atlas_elec{sbj_ix} = fn_select_elec(cfgs,atlas_elec{sbj_ix});
    roi_lab{sbj_ix}    = fn_atlas2roi_labels(atlas_elec{sbj_ix}.atlas_label,atlas_id,roi_id);
    tis_elec{sbj_ix}   = fn_select_elec(cfgs,tis_elec{sbj_ix});
    
    %% Collect SBJ results
    % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
    gm_perc{sbj_ix} = tis_elec{sbj_ix}.tissue_prob(:,1);
    gm_bin{sbj_ix}  = tis_elec{sbj_ix}.tissue_prob(:,1)>gm_thresh;
    
    % Check for ANOVA group effects + RT correlations
    grp_sig{sbj_ix} = zeros([numel(stat.label) numel(grp_lab)+1]);
    grp_sig{sbj_ix}(:,1:numel(grp_lab)) = any(qvals<0.05,3)';
    grp_sig{sbj_ix}(:,numel(grp_lab)+1) = sum(squeeze(stat.mask(:,1,:)),2)>0;
    
    clear SBJ SBJ_vars w2 stat qvals
end

%% Print positive and negative findings across analyses
% prints: ch_lab roi_lab gm% gm_bin grp1 grp2 rt
int_spacer = strjoin(repmat({'%i'},[numel(grp_lab)+2 1]),'\t');
print_line = ['%s\t\t%s\t%.2f\t' int_spacer '\n'];
grp_wm_sig_cnt = 0;
grp_sig_cnt = 0;
grp_elec_cnt = 0;
grp_wm_sig_per = [];
for sbj_ix = 1:numel(SBJs)
    wm_sig_cnt = 0;
    wm_sig_elecs = {};
    wm_sig_per = [];
    sig_cnt = 0;
    fprintf('============== %s ==============\n',SBJs{sbj_ix});
    fprintf('lab\t\troi\tgm\tgmbin\t%s\tRT\n',strjoin(grp_lab,'\t'));
    for ch_ix = 1:numel(atlas_elec{sbj_ix}.label)
        if any(grp_sig{sbj_ix}(ch_ix,:))
            sig_cnt = sig_cnt + 1;
        end
        if gm_bin{sbj_ix}(ch_ix)==0 && any(grp_sig{sbj_ix}(ch_ix,:))
            wm_sig_cnt = wm_sig_cnt + 1;
            wm_sig_elecs = [wm_sig_elecs atlas_elec{sbj_ix}.label(ch_ix)];
            wm_sig_per = [wm_sig_per gm_perc{sbj_ix}(ch_ix)];
            grp_wm_sig_cnt = grp_wm_sig_cnt + 1;
            grp_wm_sig_per = [grp_wm_sig_per gm_perc{sbj_ix}(ch_ix)];
            
            fprintf(2,print_line,atlas_elec{sbj_ix}.label{ch_ix},roi_lab{sbj_ix}{ch_ix},gm_perc{sbj_ix}(ch_ix),...
                gm_bin{sbj_ix}(ch_ix),grp_sig{sbj_ix}(ch_ix,:));
        else
%             fprintf(print_line,atlas_elec{sbj_ix}.label{ch_ix},roi_lab{sbj_ix}{ch_ix},gm_perc{sbj_ix}(ch_ix),...
%                 gm_bin{sbj_ix}(ch_ix),grp_sig{sbj_ix}(ch_ix,:));
        end
    end
    grp_elec_cnt = grp_elec_cnt + numel(atlas_elec{sbj_ix}.label);
    grp_sig_cnt = grp_sig_cnt + sig_cnt;
%     fprintf(print_line,'RT',squeeze(grp_sig(sbj_ix,:,numel(grp_lab)+1)));
%     fprintf(print_line,'Total',squeeze(gm_perc(sbj_ix,:)));
    fprintf('\nwm_sig: %i / %i elecs (%f )\n',wm_sig_cnt,numel(atlas_elec{sbj_ix}.label),...
        wm_sig_cnt/numel(atlas_elec{sbj_ix}.label));
    fprintf('\nwm_sig: %i / %i sig elecs (%f )\n',wm_sig_cnt,sig_cnt,...
        wm_sig_cnt/sig_cnt);
    fprintf('%s\n',strjoin(wm_sig_elecs,','));
    disp(wm_sig_per);
    fprintf('\n\n');
end

%% Aggregate results per ROI

%% Print elec counts per ROI across all patients
fprintf('============== GRP WM counts ==============\n');
fprintf('\ngrp_wm_sig: %i / %i elecs (%f )\n',grp_wm_sig_cnt,grp_elec_cnt,...
        grp_wm_sig_cnt/grp_elec_cnt);
disp(mean(grp_wm_sig_per));
fprintf('\ngrp_wm_sig: %i / %i sig elecs (%f )\n',grp_wm_sig_cnt,grp_sig_cnt,...
        grp_wm_sig_cnt/grp_sig_cnt);
% for roi_ix = 1:numel(roi_list)
%     fprintf('%s = %i\n',roi_list{roi_ix},sum(gm_perc(:,roi_ix)));
% end

%% Save output
if save_out
    out_dir = '/home/knight/hoycw/PRJ_Stroop/results/HFA/';    
    filename = [out_dir 'GRP_HFA_elec_cnt_ROI_' stat_id '_' an_id '.mat'];
    fprintf('Saving %s\n',filename);
    save(filename,'grp_cnt','elec_cnt');
end

end
