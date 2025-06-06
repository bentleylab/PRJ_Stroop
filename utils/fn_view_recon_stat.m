function fn_view_recon_stat(SBJ, proc_id, stat_id, an_id, view_space, reg_type, show_labels, hemi, plot_out, varargin)
%% Plot a reconstruction with electrodes colored according to statistics
%   FUTURE 1: this is a static brain, need to adapt to a movie!
%   FUTURE 2: add option for stat_var to be a cell with 2nd stat for edge
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of the stats
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'PC': ANOVA of proportion congruence (red for sig)
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

%% Process plotting params
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

% Implement the default options
ns_color = [0 0 0];
if ~exist('view_angle','var')
    if strcmp(hemi,'l')
        view_angle = [-60 30];
    elseif any(strcmp(hemi,{'r','b'}))
        view_angle = [60 30];
    else
        error(['unknown hemi: ' hemi]);
    end
end
if ~exist('mesh_alpha','var')
    if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
        mesh_alpha = 0.3;
    else
        mesh_alpha = 0.8;
    end
end
if show_labels
    lab_arg = 'label';
else
    lab_arg = 'off';
end
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

%% Load elec struct
load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_',view_space,reg_suffix,'.mat']);

%% Remove electrodes that aren't in hemisphere
if ~plot_out
    cfgs = [];
    cfgs.channel = fn_select_elec_lab_match(elec, hemi, [], []);
    elec = fn_select_elec(cfgs, elec);
end

%% Load brain recon
mesh = fn_load_recon_mesh(SBJ,view_space,reg_type,hemi);

%% Load Stats
% Determine options: {'actv','CI','RT','CNI','PC'}
if strcmp(stat_id,'actv')
    % Load HFA to get channel list
    load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '.mat']);
    % Load actv results
    load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_actv_mn100.mat']);
    elec = fn_reorder_elec(elec,hfa.label);
    grp_lab = {};
%     grp_colors = {'k','r','b'};
    elec_colors = cell([numel(hfa.label) 1]);
    for ch_ix = 1:numel(hfa.label)
        if any(strcmp(actv_ch,hfa.label{ch_ix}))
            actv_epochs = actv_ch_epochs{strcmp(actv_ch,hfa.label{ch_ix})};
            epoch_signs = zeros([size(actv_epochs,1) 1]);
            for ep_ix = 1:size(actv_epochs,1)
                sig_chunk_ix = [find(hfa.time==actv_epochs(ep_ix,1))...
                    find(hfa.time==actv_epochs(ep_ix,2))];
                % Find sign of (de)activation
                if 0<=squeeze(mean(mean(hfa.powspctrm(:,ch_ix,1,sig_chunk_ix(1):sig_chunk_ix(2)),1),4))
                    epoch_signs(ep_ix) = 1;
                else
                    epoch_signs(ep_ix) = -1;
                end
            end
            if any(epoch_signs==1) && any(epoch_signs==-1)
                elec_colors{ch_ix,1} = [0.75 0 0.75];
            elseif any(epoch_signs==1)
                elec_colors{ch_ix,1} = [1 0 0];
            elseif any(epoch_signs==-1)
                elec_colors{ch_ix,1} = [0 0 1];
            end
        else
            elec_colors{ch_ix,1} = ns_color;
        end
    end
elseif strcmp(stat_id,'CSE')
    load([SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'_',stat_id,'.mat']);
    elec = fn_reorder_elec(elec,stat.label);
    grp_lab = {};
	[~, cond_colors, ~] = fn_condition_label_styles(stat_id);
    elec_colors = cell([numel(stat.label) 1]);
    for ch_ix = 1:numel(stat.label)
        if any(stat.mask(ch_ix,1,:))
            sig_epoch = find(stat.mask(ch_ix,1,:));
            if any(diff(sig_epoch)>1)
                error([SBJ ' has multiple sig epochs in channel ' stat.label{ch_ix}]);
            end
            % Find sign of (de)activation
            cI_mean = squeeze(mean(mean(hfa{1}.powspctrm(:,ch_ix,1,sig_epoch),1),4));
            iI_mean = squeeze(mean(mean(hfa{2}.powspctrm(:,ch_ix,1,sig_epoch),1),4));
            if cI_mean>=iI_mean
                elec_colors{ch_ix,1} = cond_colors{1};
            else
                elec_colors{ch_ix,1} = cond_colors{2};
            end
        else
            elec_colors{ch_ix,1} = ns_color;
        end
    end    
else    % ANOVA
    f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
    load(f_name);
    eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
    [grp_lab, grp_colors, ~] = fn_group_label_styles(st.model_lab);
    if st.rt_corr
        [rt_lab, rt_color, ~]    = fn_group_label_styles('RT');
    end
    % % Load RTs
    % load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    
    elec = fn_reorder_elec(elec,w2.label);
    
    % FDR correct pvalues for ANOVA
%     win_lim = {}; win_center = {};
    elec_colors = cell([numel(w2.label) numel(grp_lab)+st.rt_corr]);
    for ch_ix = 1:numel(w2.label)
        % Consolidate to binary sig/non-sig
        for grp_ix = 1:numel(grp_lab)
            if any(w2.qval(grp_ix,:)<st.alpha,2)
                elec_colors{ch_ix,grp_ix} = grp_colors{grp_ix};  % sig
            else
                elec_colors{ch_ix,grp_ix} = ns_color;  % non-sig
            end
        end
        if st.rt_corr && any(stat.mask(ch_ix,1,:))
            elec_colors{ch_ix,numel(grp_lab)+st.rt_corr} = rt_color{:};  % sig
        else
            elec_colors{ch_ix,numel(grp_lab)+st.rt_corr} = ns_color;  % non-sig
        end
    end
end

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
f = {};
for grp_ix = 1:numel(grp_lab)+st.rt_corr
    if any(strcmp(stat_id,{'actv','CSE'}))
        plot_name = [SBJ '_' stat_id '_' an_id];
    else
        if grp_ix<=numel(grp_lab)
            plot_name = [SBJ '_ANOVA_' grp_lab{grp_ix} '_' stat_id '_' an_id];
        else
            plot_name = [SBJ '_ANOVA_' rt_lab{1} '_' stat_id '_' an_id];
        end
    end
    f{grp_ix} = figure('Name',plot_name);
    
    % Plot 3D mesh
    ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
    
    % Plot electrodes on top
    cfgs = [];
    for e = 1:numel(elec.label)
        cfgs.channel = elec.label(e);
        elec_tmp = fn_select_elec(cfgs, elec);
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_colors{e,grp_ix}, 'label', lab_arg);
    end
    
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
    fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
        'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
        '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
    set(f{grp_ix}, 'windowkeypressfcn',   @cb_keyboard);
end


end