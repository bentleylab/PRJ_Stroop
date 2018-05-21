function SBJ10c_HFA_GRP_elec_cnt_ROI_RT_ANOVA(SBJs,stat_id,an_id,roi_id,pipeline_id,save_out)
% Load HFA analysis results for original CI HFA test to compare to RT correlation+ANOVA results
% Original analysis:
%   'CI'- original clusterbased permutation corrected via FT results
% Compared to new analyses:
%   RT correlation: correlation with RT and cluster-based perm. corrected via FT
%   CNI ANOVA factor: FDR corrected sliding window ANOVA after regressing off RT as a confound
% OUTPUTS:
%   number of electrodes sig vs. non-sig across the two analyses,
%   separately for orig vs. RT and orig vs. CNI
% clear all; %close all;
if ischar(save_out); save_out = str2num(save_out); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Prep variables
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

[roi_list, ~, einfo_roi_col] = fn_roi_label_styles(roi_id);

[grp_lab, ~, ~] = fn_group_label_styles(model_lab);

% Set up electrode counts
grp_cnt   = zeros([numel(SBJs) numel(roi_list) numel(grp_lab)+1]);
elec_cnt  = zeros([numel(SBJs) numel(roi_list)]);

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id,'.mat'),'stat','w2');
    
    %% FDR correct pvalues for ANOVA
    qvals = NaN(size(w2.pval));
    for ch_ix = 1:numel(stat.label)
        [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
    end
    
    %% Load ROI and GM/WM info
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
    
    % Sort by gROI, then ROI
    einfo = sortrows(einfo,[3,2]);
    if ~isempty(setdiff(stat.label,einfo(:,1)))
        error('ERROR: Electrodes do not match between stat and einfo!');
    end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
        einfo_ix = strmatch(stat.label(ch_ix),einfo(:,1),'exact');
        
        % If elec matches roi_list, get stats
        if ~isempty(strmatch(einfo(einfo_ix,einfo_roi_col),roi_list,'exact'))
            roi_ix = strmatch(einfo(einfo_ix,einfo_roi_col),roi_list,'exact');
            elec_cnt(sbj_ix,roi_ix) = elec_cnt(sbj_ix,roi_ix)+1;
            
            % Check for ANOVA group effects
            for grp_ix = 1:numel(grp_lab)
                if any(squeeze(qvals(grp_ix,ch_ix,:))<0.05)
                    grp_cnt(sbj_ix,roi_ix,grp_ix) = grp_cnt(sbj_ix,roi_ix,grp_ix) + 1;
                end
            end
            
            % Check for RT correlations, get epochs
            if sum(squeeze(stat.mask(ch_ix,1,:)))>0
                grp_cnt(sbj_ix,roi_ix,numel(grp_lab)+1) = grp_cnt(sbj_ix,roi_ix,numel(grp_lab)+1) + 1;
            end
        end
    end
    clear SBJ SBJ_vars w2 stat qvals
end

%% Print positive and negative findings across analyses
int_spacer = strjoin(repmat({'%i'},size(roi_list)),'\t');
print_line = ['%s\t' int_spacer '\n'];
for sbj_ix = 1:numel(SBJs)
    fprintf('============== %s ==============\n',SBJs{sbj_ix});
    fprintf('\t%s\n',strjoin(roi_list,'\t'));
    for grp_ix = 1:numel(grp_lab)
        fprintf(print_line,grp_lab{grp_ix},squeeze(grp_cnt(sbj_ix,:,grp_ix)));
    end
    fprintf(print_line,'RT',squeeze(grp_cnt(sbj_ix,:,numel(grp_lab)+1)));
    fprintf(print_line,'Total',squeeze(elec_cnt(sbj_ix,:)));
    fprintf('\n\n');
end

%% Print elec counts per ROI across all patients
fprintf('============== GRP ROI counts ==============\n');
for roi_ix = 1:numel(roi_list)
    fprintf('%s = %i\n',roi_list{roi_ix},sum(elec_cnt(:,roi_ix)));
end

%% Save output
if save_out
    out_dir = '/home/knight/hoycw/PRJ_Stroop/results/HFA/';    
    filename = [out_dir 'GRP_HFA_elec_cnt_ROI_' stat_id '_' an_id '.mat'];
    fprintf('Saving %s\n',filename);
    save(filename,'grp_cnt','elec_cnt');
end

end
