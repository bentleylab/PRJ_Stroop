function fn_view_recon_atlas_grp_stat_venn(SBJs, proc_id, stat_conds, reg_type, show_labels,...
                                 hemi, atlas_id, roi_id, plot_out, varargin)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJs [cell array str] - subject IDs to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_conds [cell array] - {{'stat_id1','an_id1','cond1'},...,{'stat_idn','an_idn','condn'}};
%     stat_id_n [str] - ID of statistical analysis
%     cond_n [str] - ID of condition/group to overlap with
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pCNI': ANOVA of previous trial congruence (red for sig)
%       'PC': ANOVA of proportion congruence (red for sig)
%     an_id [str] - analysis ID for preprocessing, filtering, etc.
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%       'gROI','mgROI','main3' - general ROIs (lobes or broad regions)
%       'ROI','thryROI','LPFC','MPFC','OFC','INS' - specific ROIs (within these larger regions)
%       'Yeo7','Yeo17' - colored by Yeo networks
%       'tissue','tissueC' - colored by tisseu compartment, e.g., GM vs WM vs OUT
%   plot_out [0/1] - include electrodes that don't have an atlas label or in hemi?

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Handle variables
% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'view_angle')
            view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            mesh_alpha = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Define default options
% view_space = 'mni';
if ~exist('view_angle','var')
    if strcmp(hemi,'l')
        view_angle = [-90 0];
    elseif any(strcmp(hemi,{'r','b'}))
        view_angle = [90 0];
    else
        error(['unknown hemi: ' hemi]);
    end
end
if ~exist('mesh_alpha','var')
    % assume SEEG
    mesh_alpha = 0.3;
end
if show_labels
    lab_arg = 'label';
else
    lab_arg = 'off';
end
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];    % MNI space
else
    reg_suffix = '';                % Patient space
end

if numel(stat_conds) < 2 || numel(stat_conds) > 3; error('why venn?'); end
for st_ix = 2:numel(stat_conds)
    if numel(stat_conds{st_ix})~=3; error('need stat_id, an_id, cond in each stat_cond'); end
end

% Venn colors
venn_colors = fn_venn_colors(numel(stat_conds));
all_color   = [1 1 1];

%% Prep report
% Organize IDs
stat_ids = cell(size(stat_conds)); an_ids = cell(size(stat_conds)); cond_ids = cell(size(stat_conds));
for st_ix = 1:numel(stat_conds)
    stat_ids{st_ix} = stat_conds{st_ix}{1};
    an_ids{st_ix}   = stat_conds{st_ix}{2};
    cond_ids{st_ix} = stat_conds{st_ix}{3};
end
% Create output dir
save_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_overlap/' strjoin(stat_ids,'-')...
            '/' strjoin(an_ids,'-') '/' strjoin(cond_ids,'-') '/'];
if ~exist(save_dir,'dir')
    [~] = mkdir(save_dir);
end
% Create report
sig_report_fname = [save_dir 'GRP_' atlas_id '_' roi_id '_' hemi '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');
% Print IDs
title_str = repmat('%-40s',[1 numel(stat_conds)]);
fprintf(sig_report,[title_str '\n'],stat_ids{:});
fprintf(sig_report,[title_str '\n'],an_ids{:});
fprintf(sig_report,[title_str '\n'],cond_ids{:});
fprintf(sig_report,[repmat('=',[1 40*numel(stat_conds)]) '\n']);
% Print result headings
fprintf(sig_report,[repmat('%-20s',[1 numel(stat_conds)+2]) '\n'],'label','ROI',cond_ids{:});
fprintf(sig_report,[repmat('-',[1 40*numel(stat_conds)]) '\n']);
result_str = ['%-20s%-20s' repmat('%-20i',[1 numel(stat_conds)]) '\n'];

%% Load Data
elec_sbj = cell([numel(SBJs) 1]);
good_sbj = true([numel(SBJs) 1]);
sig_mat  = cell([numel(SBJs) 1]);
% all_roi_labels = cell([numel(stat_conds) 1]);
all_roi_colors = cell([numel(SBJs) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %% Prepare elec structs
    % Load elec struct
    try
        elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',reg_suffix,'_',atlas_id,'_full.mat'];
        if exist([elec_fname(1:end-4) '_' roi_id '.mat'],'file')
            elec_fname = [elec_fname(1:end-4) '_' roi_id '.mat'];
        end
        tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
    catch
        error([elec_fname 'doesnt exist, exiting...']);
    end
    
    % Append SBJ name to labels
    orig_labels = elec_sbj{sbj_ix}.label;
    for e_ix = 1:numel(elec_sbj{sbj_ix}.label)
        elec_sbj{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix}.label{e_ix}];
    end
    
    % Match elecs to atlas ROIs
    if any(strcmp(atlas_id,{'DK','Dx','Yeo7'}))
        if ~isfield(elec_sbj{sbj_ix},'man_adj')
            elec_sbj{sbj_ix}.roi       = fn_atlas2roi_labels(elec_sbj{sbj_ix}.atlas_lab,atlas_id,roi_id);
        end
        if strcmp(roi_id,'tissueC')
            elec_sbj{sbj_ix}.roi_color = fn_tissue2color(elec_sbj{sbj_ix});
        elseif strcmp(atlas_id,'Yeo7')
            elec_sbj{sbj_ix}.roi_color = fn_atlas2color(atlas_id,elec_sbj{sbj_ix}.roi);
        else
            elec_sbj{sbj_ix}.roi_color = fn_roi2color(elec_sbj{sbj_ix}.roi);
        end
    elseif any(strcmp(atlas_id,{'Yeo17'}))
        if ~isfield(elec_sbj{sbj_ix},'man_adj')
            elec_sbj{sbj_ix}.roi       = elec_sbj{sbj_ix}.atlas_lab;
        end
        elec_sbj{sbj_ix}.roi_color = fn_atlas2color(atlas_id,elec_sbj{sbj_ix}.roi);
    end
    
    % Remove hemi and/or atlas elecs
    if ~plot_out
        % Remove electrodes that aren't in atlas ROIs & hemisphere
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, atlas_id, roi_id);
    else
        % Remove electrodes that aren't in hemisphere
        roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, [], []);
    end
    cfgs = []; cfgs.channel = roi_elecs;
    elec_sbj{sbj_ix} = fn_select_elec(cfgs,elec_sbj{sbj_ix});
    
    %% Load Stats
    if ~isempty(roi_elecs)
        sig_mat{sbj_ix} = zeros([numel(elec_sbj{sbj_ix}.label) numel(stat_conds)]);
        for st_ix = 1:numel(stat_conds)
            if strcmp(cond_ids{st_ix},'actv')
                error('not done actv yet');
                % load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_' stat_id '.mat'],'actv');
            elseif strcmp(cond_ids{st_ix},'RT')
                error('RT not implemented');
            else    % ANOVA
                load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_ids{st_ix} '_' an_ids{st_ix} '.mat']);
                for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
                    stat_ch_ix = strcmp(w2.label,orig_labels{ch_ix});
                    if ~any(stat_ch_ix)
                        fprintf(2,'\tWARNING: elec chan %s not found in w2!\n',elec_sbj{sbj_ix}.label{ch_ix});
                    end
                    % Consolidate to binary sig/non-sig
                    if any(squeeze(w2.qval(strcmp(st.groups,cond_ids{st_ix}),stat_ch_ix,:))<st.alpha)
                        sig_mat{sbj_ix}(ch_ix,st_ix) = 1;
                    end
                end
                clear w2 st actv
            end
        end
        
        % Report on significant electrodes for this SBJ (even if non-sig)
        for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
            fprintf(sig_report,result_str,elec_sbj{sbj_ix}.label{ch_ix},...
                                          elec_sbj{sbj_ix}.roi{ch_ix},sig_mat{sbj_ix}(ch_ix,:));
        end
        
        %% Compile Statistics
        if any(sig_mat{sbj_ix}(:))
            print_nums = zeros([2 numel(stat_conds)]);
            for st_ix = 1:numel(stat_conds)
                print_nums(1,st_ix) = sum(sig_mat{sbj_ix}(:,st_ix));
                print_nums(2,st_ix) = sum(sig_mat{sbj_ix}(:,st_ix))/size(sig_mat{sbj_ix},1);
            end
            fprintf(['%-10s' repmat('%-3i(%.3f)\t',[1 numel(stat_conds)]) '\n'],[SBJ ' sig:'],print_nums(:));
            
            % Select color and add to elec
            all_roi_colors{sbj_ix} = zeros([sum(any(sig_mat{sbj_ix},2)) 3]);
            sig_elecs = cell([sum(any(sig_mat{sbj_ix},2)) 1]);
            sig_ix = 0;
            for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
                if any(sig_mat{sbj_ix}(ch_ix,:))
                    sig_ix = sig_ix + 1;
                    sig_elecs{sig_ix} = elec_sbj{sbj_ix}.label{ch_ix};
                    % Select color
                    cond_ix = find(sig_mat{sbj_ix}(ch_ix,:));
                    if numel(cond_ix)==numel(stat_conds) && numel(stat_conds) == 3
                        all_roi_colors{sbj_ix}(sig_ix,:) = all_color;
                    elseif numel(cond_ix)==2
                        all_roi_colors{sbj_ix}(sig_ix,:) = venn_colors{cond_ix(1),cond_ix(2)};
                    else
                        all_roi_colors{sbj_ix}(sig_ix,:) = venn_colors{cond_ix,cond_ix};
                    end
                end
            end
            
            % Select sig elecs for plotting
            cfgs = []; cfgs.channel = sig_elecs;
            elec_sbj{sbj_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix});
        else
            % Print no significant elecs
            elec_sbj{sbj_ix} = {}; good_sbj(sbj_ix) = false;
            fprintf(2,'\t%s has no significant channels for any stat\n',SBJ);
        end
        
    else
        % Print no ROI match
        elec_sbj{sbj_ix} = {}; good_sbj(sbj_ix) = false;
        fprintf(2,'\t%s has no channels in %s hemi %s\n',SBJ,atlas_id,hemi);
    end
    clear SBJ SBJ_vars SBJ_vars_cmd sig_elecs sig_ix cond_ix 
end

%% Finish report and save sig_mat
fclose(sig_report);
sig_elecs_fname = [save_dir 'GRP_' atlas_id '_' roi_id '_sig_elecs.mat'];
save(sig_elecs_fname,'-v7.3','stat_conds','SBJs','elec_sbj','sig_mat');

%% Combine elec structs
elec = ft_appendsens([],elec_sbj{good_sbj});
elec.roi_color = vertcat(all_roi_colors{:});    % appendsens strips that field
%     elec.roi       = all_roi_labels{cond_ix};    % appendsens strips that field

%% Load brain recon
mesh = fn_load_recon_mesh([],'mni',reg_type,hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
% Plot venn recon
plot_name = ['GRP_' strjoin(stat_ids,'-')];
f = figure('Name',plot_name);

% Plot 3D mesh
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);

% Plot electrodes on top
for e = 1:numel(elec.label)
    cfgs = []; cfgs.channel = elec.label(e);
    elec_tmp = fn_select_elec(cfgs,elec);
    ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.roi_color, 'label', lab_arg);
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(f, 'windowkeypressfcn',   @cb_keyboard);

%% Plot legend in separate figure
figure; hold on;
venn_legend = {}; leg_ix = 0;
for cond_ix1 = 1:numel(stat_conds);
    for cond_ix2 = cond_ix1:numel(stat_conds);
        scatter(cond_ix1,cond_ix2,500,'MarkerFaceColor',venn_colors{cond_ix1,cond_ix2},...
                'MarkerEdgeColor','k');
        leg_ix = leg_ix + 1;
        if cond_ix1==cond_ix2
            venn_legend{leg_ix} = [cond_ids{cond_ix1} num2str(cond_ix1)];
        else
            venn_legend{leg_ix} = [cond_ids{cond_ix1} num2str(cond_ix1) '+'...
                                   cond_ids{cond_ix2} num2str(cond_ix2)];
        end
    end
end
if numel(stat_conds)>2
    scatter(numel(stat_conds),1,500,'MarkerFaceColor','w','MarkerEdgeColor','k');
    venn_legend{leg_ix+1} = strjoin(cond_ids,'+');
end
xlim([0 numel(stat_conds)+1]);
ylim([0 numel(stat_conds)+1]);
legend(venn_legend);

end
