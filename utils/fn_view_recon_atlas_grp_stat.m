function fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_id, reg_type, show_labels,...
                                 hemi, atlas_id, roi_id, mirror, plot_out, save_fig, fig_ftype, varargin)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJs [cell array str] - subject IDs to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of the stats
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pCNI': ANOVA of previous trial congruence (red for sig)
%       'PC': ANOVA of proportion congruence (red for sig)
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
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
% Add default view_angle if not defined
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

% ROI info
[roi_list, ~] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end
fprintf('Using atlas: %s\n',atlas_id);

%% Prep report
out_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_recons/' stat_id '/' an_id '/'];
if ~exist(out_dir,'dir')
    [~,~] = mkdir(out_dir);
end
sig_report_fname = [out_dir 'GRP_' atlas_id '_' roi_id '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');

%% Load Data
% Get condition info
load([root_dir,'PRJ_Stroop/data/',SBJs{1},'/04_proc/',SBJs{1},'_smANOVA_ROI_',stat_id,'_',an_id,'.mat'],'st');
cond_lab = st.groups;

elec_sbj = cell([numel(SBJs) numel(cond_lab)]);
good_sbj = true([numel(SBJs) numel(cond_lab)]);
all_roi_labels = cell([numel(cond_lab) 1]);
all_roi_colors = cell([numel(cond_lab) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('\tLoading %s (%i / %i)...\n',SBJ,sbj_ix,numel(SBJs));
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load elec struct
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',reg_suffix,'_',atlas_id,'_final.mat'];
    tmp = load(elec_fname); elec_sbj{sbj_ix,1} = tmp.elec;
    
    % Append SBJ name to labels
    elec_sbj{sbj_ix,1}.color = cell(size(elec_sbj{sbj_ix,1}.label));
    for e_ix = 1:numel(elec_sbj{sbj_ix,1}.label)
        elec_sbj{sbj_ix,1}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix,1}.label{e_ix}];
        elec_sbj{sbj_ix,1}.color{e_ix} = fn_roi2color(elec_sbj{sbj_ix,1}.(roi_field){e_ix});
    end
        
    % Remove hemi and/or atlas elecs
    if ~plot_out
        % Remove electrodes that aren't in atlas ROIs & hemisphere
        if mirror
            roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix,1}, 'b', atlas_id, roi_id);
        else
            roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix,1}, hemi, atlas_id, roi_id);
        end
    else
        % Remove electrodes that aren't in hemisphere
        if mirror
            roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix,1}, 'b', [], []);
        else
            roi_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix,1}, hemi, [], []);
        end
    end
    
    % Mirror hemispheres
    if mirror
        elec_sbj{sbj_ix,1}.chanpos(~strcmp(elec_sbj{sbj_ix,1}.hemi,hemi),1) = ...
                        -elec_sbj{sbj_ix,1}.chanpos(~strcmp(elec_sbj{sbj_ix,1}.hemi,hemi),1);
        hemi_str = [hemi 'b'];
    else
        hemi_str = hemi;
    end
    
    % Copy for other conditions
    for cond_ix = 2:numel(cond_lab)
        elec_sbj{sbj_ix,cond_ix} = elec_sbj{sbj_ix,1};
    end
    
    % Load Stats
    % Determine options: {'actv','CI','RT','CNI','PC'}
    sig_ch = cell(size(cond_lab));
    if contains(stat_id,'actv')
        load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_' stat_id '.mat'],'actv');
        sig_ch{1} = actv.label(any(actv.mask,2));
        clear actv
    else    % ANOVA
        load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id '_' an_id '.mat'],'w2');
        
        for ch_ix = 1:numel(w2.label)
            % Consolidate to binary sig/non-sig
            for cond_ix = 1:numel(cond_lab)
                if any(squeeze(w2.qval(cond_ix,ch_ix,:))<=st.alpha)
                    sig_ch{cond_ix} = [sig_ch{cond_ix} {[SBJs{sbj_ix} '_' w2.label{ch_ix}]}];
                end
            end
        end
        clear w2
    end
    
    % Select sig elecs
    fprintf(sig_report,'===============================================================\n');
    for cond_ix = 1:numel(cond_lab)
        % Report on significant electrodes for this SBJ
        fprintf(sig_report,'\t%s - %s  = %i / %i (%.02f) sig elecs:\n',SBJs{sbj_ix},cond_lab{cond_ix},...
            numel(sig_ch{cond_ix}),numel(elec_sbj{sbj_ix,cond_ix}.label),...
            100*numel(sig_ch{cond_ix})/numel(elec_sbj{sbj_ix,cond_ix}.label));
        for sig_ix = 1:numel(sig_ch{cond_ix})
            e_ix = find(strcmp(elec_sbj{sbj_ix,cond_ix}.label,sig_ch{cond_ix}{sig_ix}));
            fprintf(sig_report,'%s - %s (%s)\n',sig_ch{cond_ix}{sig_ix},elec_sbj{sbj_ix,cond_ix}.(roi_field){e_ix},...
                elec_sbj{sbj_ix,cond_ix}.hemi{e_ix});
        end
        
        % Select sig elecs && elecs matching atlas
        % fn_select_elec messes up if you try to toss all elecs
        good_elecs = intersect(roi_elecs, sig_ch{cond_ix});
        if numel(intersect(elec_sbj{sbj_ix,cond_ix}.label,good_elecs))==0
            elec_sbj{sbj_ix,cond_ix} = {};
            good_sbj(sbj_ix,cond_ix) = false;
            warning('WARNING!!! All sig_ch are out of hemisphere and/or ROI!');
            fprintf(sig_report,'\t!!! 0 sig_ch remain\n');
        else
            cfgs = [];
            cfgs.channel = good_elecs;
            elec_sbj{sbj_ix,cond_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix,cond_ix});
            all_roi_labels{cond_ix} = [all_roi_labels{cond_ix}; elec_sbj{sbj_ix,cond_ix}.(roi_field)];
            all_roi_colors{cond_ix} = [all_roi_colors{cond_ix}; elec_sbj{sbj_ix,cond_ix}.color];
            fprintf(sig_report,'\t%i sig_ch remain\n',numel(elec_sbj{sbj_ix,cond_ix}.label));
        end
        fprintf(sig_report,'---------------------------------------------------------------\n');
    end
    fprintf(sig_report,'===============================================================\n');
    clear SBJ SBJ_vars SBJ_vars_cmd
end

%% Combine elec structs
elec = cell([numel(cond_lab) 1]);
for cond_ix = 1:numel(cond_lab)
    elec{cond_ix} = ft_appendsens([],elec_sbj{good_sbj(:,cond_ix),cond_ix});
    elec{cond_ix}.(roi_field) = all_roi_labels{cond_ix};    % appendsens strips that field
    elec{cond_ix}.color       = all_roi_colors{cond_ix};    % appendsens strips that field
end

%% Load brain recon
mesh = fn_load_recon_mesh([],'mni',reg_type,'pial',hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
f = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    if contains(stat_id,'actv')
        plot_name = ['GRP_' stat_id '_' an_id '_' hemi_str '_' view_str];
    else
        plot_name = ['GRP_' cond_lab{cond_ix} '_' stat_id '_' an_id '_' hemi_str '_' view_str];
    end
    f{cond_ix} = figure('Name',plot_name);
    
    % Plot 3D mesh
    ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
    
    % Plot electrodes on top
    for e = 1:numel(elec{cond_ix}.label)
        cfgs = []; cfgs.channel = elec{cond_ix}.label(e);
        elec_tmp = fn_select_elec(cfgs,elec{cond_ix});
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.color{1}, 'label', lab_arg);
    end
    
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
    fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
        'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
        '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
    set(f{cond_ix}, 'windowkeypressfcn',   @cb_keyboard);
    
    if save_fig
        saveas(f{cond_ix}, [out_dir plot_name '.' fig_ftype]);
    end
end

end
