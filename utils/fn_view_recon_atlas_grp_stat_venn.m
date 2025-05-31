function fn_view_recon_atlas_grp_stat_venn(SBJs, proc_id, stat_conds, reg_type, show_labels,...
                                 hemi, atlas_id, roi_id, mirror, plot_out, save_fig, fig_ftype, varargin)
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
%   mirror [0/1] - plot the other hemi, 
%   plot_out [0/1] - include electrodes that don't have an atlas label or in hemi?
%   save_fig [0/1] - save the .fig file?
%   fig_ftype [str] - extension of saved figure

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
    view_str = 'def';
end
% Adjust view angle if custom
if ischar(view_angle)
    view_str = view_angle;
    if (strcmp(hemi,'l') && strcmp(view_angle,'med')) || (strcmp(hemi,'r') && strcmp(view_angle,'lat'))
        view_angle = [90 0];
    elseif (strcmp(hemi,'l') && strcmp(view_angle,'lat')) || (strcmp(hemi,'r') && strcmp(view_angle,'med'))
        view_angle = [-90 0];
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
[roi_list, ~] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end

%% Organize IDs
if numel(stat_conds) < 2 || numel(stat_conds) > 3; error('why venn?'); end
for st_ix = 2:numel(stat_conds)
    if numel(stat_conds{st_ix})~=3; error('need stat_id, an_id, cond in each stat_cond'); end
end

stat_ids = cell(size(stat_conds)); an_ids = cell(size(stat_conds)); cond_ids = cell(size(stat_conds));
for st_ix = 1:numel(stat_conds)
    stat_ids{st_ix} = stat_conds{st_ix}{1};
    an_ids{st_ix}   = stat_conds{st_ix}{2};
    cond_ids{st_ix} = stat_conds{st_ix}{3};
end

% Venn colors
if strcmp(an_ids{1},'HGm_S2t151_zbtA_sm0_wn100') && ...
        strcmp(an_ids{2},'HGm_S2t251_zbtA_sm0_wn100') && ...
        strcmp(an_ids{3},'HGm_R5t101_zbtA_sm0_wn100') && ...
        strcmp(cond_ids{1},cond_ids{2}) && strcmp(cond_ids{2},cond_ids{3})
    color_id = 'SDR';
elseif strcmp(cond_ids{1},'CNI') && strcmp(cond_ids{2},'pCNI') && strcmp(cond_ids{3},'PC')
    color_id = 'ConCSPC';
elseif strcmp(cond_ids{1},'pCNI') && strcmp(cond_ids{2},'PC')
    color_id = 'CSPC';
end
if exist('color_id','var')
    venn_colors = fn_venn_colors(numel(stat_conds), 'cond_id', color_id);
else
    venn_colors = fn_venn_colors(numel(stat_conds));
end
all_color   = [1 1 1];

%% Load Data
elec_sbj    = cell([numel(SBJs) 1]);
elec_sig    = cell([numel(SBJs) 1]);
good_sbj    = true([numel(SBJs) 1]);
sig_roi_mat = cell([numel(SBJs) 1]);
roi_mat     = cell([numel(SBJs) 1]);
all_roi_colors = cell([numel(SBJs) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %% Prepare elec structs
    % Load elec struct
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',reg_suffix,'_',atlas_id,'_final.mat'];
    tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
    
    % Append SBJ name to labels
    orig_labels = elec_sbj{sbj_ix}.label;
    for e_ix = 1:numel(elec_sbj{sbj_ix,1}.label)
        elec_sbj{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix}.label{e_ix}];
    end
    
    % Remove hemi and/or atlas elecs
    if ~plot_out
        % Remove electrodes that aren't in atlas ROIs & hemisphere
        if mirror
            roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, 'b', atlas_id, roi_id);
            hemi_str = [hemi 'b'];
        else
            roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, atlas_id, roi_id);
            hemi_str = hemi;
        end
    else
        % Remove electrodes that aren't in hemisphere
        if mirror
            roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, 'b', [], []);
            hemi_str = [hemi 'b'];
        else
            roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, [], []);
            hemi_str = hemi;
        end
    end
    
    roi_mat{sbj_ix} = zeros([numel(elec_sbj{sbj_ix}.label) 1]);
    for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
        if any(strcmp(elec_sbj{sbj_ix}.label{ch_ix},roi_elecs))
            roi_mat{sbj_ix}(ch_ix) = find(strcmp(roi_list,elec_sbj{sbj_ix}.(roi_field){ch_ix}));
        end
    end
    
    %% Load Stats
    sig_mat = zeros([numel(elec_sbj{sbj_ix}.label) numel(stat_conds)]);
    if any(roi_mat{sbj_ix})
        % Mirror hemispheres
        if mirror
            elec_sbj{sbj_ix}.chanpos(~strcmp(elec_sbj{sbj_ix}.hemi,hemi),1) = ...
                -elec_sbj{sbj_ix}.chanpos(~strcmp(elec_sbj{sbj_ix}.hemi,hemi),1);
        end
        
        for st_ix = 1:numel(stat_conds)
            % Load correct stat output structure, ensure match to elec
            if strcmp(cond_ids{st_ix},'actv')
                error('not done actv yet');
                % load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_' stat_id '.mat'],'actv');
            elseif strcmp(cond_ids{st_ix},'RT')
                error('RT not implemented');
            else    % ANOVA
                load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_ids{st_ix} '_' an_ids{st_ix} '.mat']);
                if numel(elec_sbj{sbj_ix}.label)~=numel(w2.label) || ~all(strcmp(orig_labels,w2.label))
                    error('mismatched labels in elec and w2!');
                end
            end
            % Consolidate to binary sig/non-sig
            for ch_ix = 1:numel(elec_sbj{sbj_ix}.label)
                if any(squeeze(w2.qval(strcmp(st.groups,cond_ids{st_ix}),ch_ix,:))<=st.alpha)
                    sig_mat(ch_ix,st_ix) = 1;
                end
            end
            clear w2 st rt actv
        end
        
        %% Compile Statistics
        if any(any(sig_mat(roi_mat{sbj_ix}~=0,:),2))
            % Print # and % sig
            print_nums = zeros([2 numel(stat_conds)]);
            for st_ix = 1:numel(stat_conds)
                print_nums(1,st_ix) = sum(sig_mat(:,st_ix));
                print_nums(2,st_ix) = sum(sig_mat(:,st_ix))/size(sig_mat,1);
            end
            fprintf(['%-10s' repmat('%-3i(%.3f)\t',[1 numel(stat_conds)]) '\n'],[SBJ ' sig:'],print_nums(:));
            
            % Keep intersection of significant and ROI matched electrodes
            sig_idx = any(sig_mat,2);
            roi_idx = roi_mat{sbj_ix}~=0;
            sig_roi_mat{sbj_ix} = sig_mat(roi_idx & sig_idx,:).*roi_mat{sbj_ix}(roi_idx & sig_idx);
            % Select sig elecs for plotting
            cfgs = []; cfgs.channel = elec_sbj{sbj_ix}.label(sig_idx & roi_idx);
            elec_sig{sbj_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix});
            
            % Find colors for sig elecs
            all_roi_colors{sbj_ix} = zeros([numel(elec_sig{sbj_ix}.label) 3]);
            for sig_ix = 1:numel(elec_sig{sbj_ix}.label)
                cond_ix = find(sig_roi_mat{sbj_ix}(sig_ix,:));
                if numel(cond_ix)==numel(stat_conds) && numel(stat_conds) == 3
                    all_roi_colors{sbj_ix}(sig_ix,:) = all_color;
                elseif numel(cond_ix)==2
                    all_roi_colors{sbj_ix}(sig_ix,:) = venn_colors{cond_ix(1),cond_ix(2)};
                else
                    all_roi_colors{sbj_ix}(sig_ix,:) = venn_colors{cond_ix,cond_ix};
                end
            end
            
            fprintf('\t%s has %i sig channels in %s hemi %s\n',SBJ,size(sig_roi_mat{sbj_ix},1),atlas_id,hemi);
        else
            good_sbj(sbj_ix) = false;
            fprintf(2,'\t%s has %i channels in %s hemi %s, but none are significant\n',...
                SBJ,sum(roi_mat{sbj_ix}~=0),atlas_id,hemi);
        end
    else
        % Print no ROI match
        good_sbj(sbj_ix) = false;
        fprintf(2,'\t%s has no channels in %s hemi %s\n',SBJ,atlas_id,hemi_str);
    end
    clear SBJ SBJ_vars SBJ_vars_cmd sig_ix cond_ix sig_mat roi_idx sig_idx
end

%% Combine elec structs
elec = ft_appendsens([],elec_sig{good_sbj});
elec.color = vertcat(all_roi_colors{:});    % appendsens strips that field
%     elec.roi       = all_roi_labels{cond_ix};    % appendsens strips that field

%% Load brain recon
mesh = fn_load_recon_mesh([],'mni',reg_type,'pial',hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
out_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_venn_recons/' strjoin(stat_ids,'-') '/'...
           strjoin(an_ids,'-') '/' strjoin(cond_ids,'-') '/'];
if ~exist(out_dir,'dir')
    [~,~] = mkdir(out_dir);
end
% Plot venn recon
plot_name = ['GRP_' strjoin(stat_ids,'-') '_' hemi_str '_' view_str];
f = figure('Name',plot_name);

% Plot 3D mesh
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);

% Plot electrodes on top
for e = 1:numel(elec.label)
    cfgs = []; cfgs.channel = elec.label(e);
    elec_tmp = fn_select_elec(cfgs,elec);
    ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.color, 'label', lab_arg);
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(f, 'windowkeypressfcn',   @cb_keyboard);

if save_fig
    saveas(f, [out_dir plot_name '.' fig_ftype]);
end

%% Plot legend in separate figure
venn_name = [plot_name '_leg'];
v = figure('Name',venn_name); hold on;
venn_legend = {}; leg_ix = 0;
for cond_ix1 = 1:numel(stat_conds)
    for cond_ix2 = cond_ix1:numel(stat_conds)
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

if save_fig
    saveas(v, [out_dir venn_name '.' fig_ftype]);
end

end
