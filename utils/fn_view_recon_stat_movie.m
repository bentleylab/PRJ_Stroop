function fn_view_recon_stat_movie(SBJ, proc_id, stat_id, an_id, view_space, reg_type, hemi, plot_out, plt_id, varargin)
error('fix for new st.ANOVA params');
%% Plot a reconstruction with electrodes colored according to statistics
%   FUTURE 2: add option for stat_var to be a cell with 2nd stat for edge
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of the stats
%     current:
%       'actv': HFA activation vs. baseline; red for active, blue for deactive, yellow for both
%     future:
%       'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'PC': ANOVA of proportion congruence (red for sig)
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%       '' (or anything not v/s) for patient space
%   hemi [str] - {'r','l','b'} for right, left, or both hemispheres (shouldn't be both here, too hard to see...)
%   plot_out [0/1] - plot electrodes outside of the hemisphere or ROIs of interest?
%   plt_id [str] - which set of plotting params
%   varargin [cell]:
%       view_angle [int,int] - angle to view the brain from
%       mesh_alpha [float] - transparency of the surface mesh

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

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
if ~exist('view_angle','var')
    if strcmp(hemi,'r')
        view_angle = [60 30];
    elseif strcmp(hemi,'l')
        view_angle = [-60 30];
    else
        error(['unknown hemi: ' hemi]);
    end
elseif strcmp(hemi,'b')
    error('no sense in making a both hemi movie!');
end
if ~exist('mesh_alpha','var')
    if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
        mesh_alpha = 0.3;
    else
        mesh_alpha = 0.8;
    end
end
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

% SHEILA: move these to plt_vars
ns_color = [0 0 0];
vid_ext = '.mp4';
vid_encoding = 'MPEG-4';

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
    hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
    load(hfa_fname);

    load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_actv_mn100.mat']);
    cfgs = []; cfgs.channel = actv_ch;
    hfa = ft_selectdata(cfgs,hfa);
    elec = fn_select_elec(cfgs,elec);
    
    plot_dat = squeeze(mean(hfa.powspctrm,1));
    % get clim
    clim = [prctile(abs(hfa.powspctrm(:)),5) prctile(abs(hfa.powspctrm(:)),95)];
    
    sig_t_ix = cell(size(hfa.label));
    for e = 1:numel(hfa.label)
        for ep = 1:size(actv_ch_epochs{e},1)
            sig_t_ix{e} = [sig_t_ix{e} find(hfa.time==actv_ch_epochs{e}(ep,1)):find(hfa.time==actv_ch_epochs{e}(ep,2))];
        end
    end
% elseif strcmp(stat_id,'CSE')
%     error('no CSE yet');
%     load([SBJ_vars.dirs.proc,SBJ,'_',stat_id,'_ROI_',an_id,'.mat']);
%     elec = fn_reorder_elec(elec,stat.label);
%     grp_lab = {};
% 	[~, cond_colors, ~] = fn_condition_label_styles(stat_id);
%     elec_colors = cell([numel(stat.label) 1]);
%     for ch_ix = 1:numel(stat.label)
%         if any(stat.mask(ch_ix,1,:))
%             sig_epoch = find(stat.mask(ch_ix,1,:));
%             if any(diff(sig_epoch)>1)
%                 error([SBJ ' has multiple sig epochs in channel ' stat.label{ch_ix}]);
%             end
%             % Find sign of (de)activation
%             cI_mean = squeeze(mean(mean(hfa{1}.powspctrm(:,ch_ix,1,sig_epoch),1),4));
%             iI_mean = squeeze(mean(mean(hfa{2}.powspctrm(:,ch_ix,1,sig_epoch),1),4));
%             if cI_mean>=iI_mean
%                 elec_colors{ch_ix,1} = cond_colors{1};
%             else
%                 elec_colors{ch_ix,1} = cond_colors{2};
%             end
%         else
%             elec_colors{ch_ix,1} = ns_color;
%         end
%     end    
% else    % ANOVA
%     error('no ANOVA yet)');
%     eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
%     [grp_lab, grp_colors, ~] = fn_group_label_styles(model_lab);
%     [rt_lab, rt_color, ~]    = fn_group_label_styles('RT');
%     % % Load RTs
%     % load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
%     
%     f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
%     load(f_name,'stat','w2');
%     elec = fn_reorder_elec(elec,stat.label);
%     
%     % FDR correct pvalues for ANOVA
% %     win_lim = {}; win_center = {};
%     elec_colors = cell([numel(w2.label) numel(grp_lab)+1]);
%     for ch_ix = 1:numel(stat.label)
%         pvals = squeeze(w2.pval(:,ch_ix,:));
%         [~, ~, ~, qvals] = fdr_bh(pvals);%,0.05,'pdep','yes');
%         
%         % Consolidate to binary sig/non-sig
%         for grp_ix = 1:numel(grp_lab)
%             if any(qvals(grp_ix,:)<0.05,2)
%                 elec_colors{ch_ix,grp_ix} = grp_colors{grp_ix};  % sig
%             else
%                 elec_colors{ch_ix,grp_ix} = ns_color;  % non-sig
%             end
%         end
%         if any(stat.mask(ch_ix,1,:))
%             elec_colors{ch_ix,numel(grp_lab)+1} = rt_color{:};  % sig
%         else
%             elec_colors{ch_ix,numel(grp_lab)+1} = ns_color;  % non-sig
%         end
%     end
%     
% %     % Get Sliding Window Parameters
% %     win_lim    = fn_sliding_window_lim(stat.time,win_len,win_step);
% %     win_center = round(mean(win_lim,2));
%     
% %     % Convert % explained variance to 0-100 scale
% %     w2.trial = w2.trial*100;
%     
%     
% %     % Trim data to plotting epoch
% %     %   NOTE: stat should be on stat_lim(1):stat_lim(2)+0.001 time axis
% %     %   w2 should fit within that since it's averaging into a smaller window
% %     cfg_trim = [];
% %     cfg_trim.latency = plt_vars.plt_lim_SR;
% %     stat = ft_selectdata(cfg_trim,stat);
end

%% Load timing info
evnt_times = cell(size(plt_vars.evnt_type));
for evnt_ix = 1:numel(plt_vars.evnt_type)
    if strcmp(plt_vars.evnt_type{evnt_ix},'S')
        evnt_times{evnt_ix} = 0:1/sample_rate:plt_vars.evnt_len;
    elseif strcmp(plt_vars.evnt_lab{evnt_ix},'R')
        load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
        evnt_times{evnt_ix} = mean(trial_info.response_time):1/sample_rate:mean(trial_info.response_time)+plt_vars.evnt_len;
        %     % Compile cond_type, RT, trial_n
        %     [cond_lab, cond_colors, ~] = fn_condition_label_styles('CI');
        %     cond_RTs = zeros(size(cond_lab));
        %     for cond_ix = 1:length(cond_lab)
        %         cond_RTs(cond_ix) = mean(trial_info.response_time(logical(fn_condition_index(cond_lab{cond_ix},...
        %             trial_info.condition_n))));
        %     end
        %     n_
    end
end

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
if any(strcmp(stat_id,{'actv','CSE'}))
    plot_name = [SBJ '_' stat_id '_' an_id];
else
    if grp_ix<=numel(grp_lab)
        plot_name = [SBJ '_ANOVA_' grp_lab{grp_ix} '_' stat_id '_' an_id];
    else
        plot_name = [SBJ '_ANOVA_' rt_lab{1} '_' stat_id '_' an_id];
    end
end
f = figure('Name',plot_name); hold on;

% Plot 3D mesh
mesh_obj = ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
%SHEILA: this would be the first of potentially more than one view_angle
view(view_angle); material dull; lighting gouraud;

[xsp, ysp, zsp] = sphere(100);


% Plot Timing info
% Text for time in sec (X,Y,Z in mm)
time_str = text(plt_vars.time_str_pos(1), plt_vars.time_str_pos(2), plt_vars.time_str_pos(3),...
    [num2str(hfa.time(1),'%.2f') ' s'], 'color', plt_vars.time_str_color,...
    'fontsize', plt_vars.time_str_size, 'horizontalalignment', plt_vars.time_str_horzalign);
evnt_str = text(plt_vars.evnt_str_pos(1), plt_vars.evnt_str_pos(2), plt_vars.evnt_str_pos(3),...
    plt_vars.evnt_lab{1}, 'color', plt_vars.evnt_str_color,...
    'fontsize', plt_vars.evnt_str_size, 'horizontalalignment', plt_vars.evnt_str_horzalign,...
    'fontweight', plt_vars.evnt_str_weight);
set(evnt_str,'Visible','off');

%% Plot frames
frames = struct('cdata',[],'colormap',[]);%cell(size(hfa.time));
for t_ix = 1:numel(hfa.time)
    % Plot electrodes on top
    for e = 1:numel(elec.label)
        % Sphere settings
        %!!! SHEILA: move to plt_vars
        reg_sz = 0.5;
        cmap = colormap('parula');
        elec_cmap = linspace(clim(1),clim(2),size(cmap,1));
        
        sz_lim = [0.5 5];       % radius (in mm) of spheres
        % SHEILA STARTS HERE: sz_map = linspace(sz_lim(1),sz_lim(2), size(cmap,1));
        e_sphr_ns = cell(size(elec.label));
        % e_sphr_s    = cell(size(elec.label));
        for e = 1:numel(elec.label)
            % Create sphere
            if any(sig_times{e}==t_ix)
%                 <pick the right size>
            else
%                 <non-sig size>
            end
            e_sphr_ns{e} = surf(sz_lim(1)*xsp+elec.chanpos(e,1), sz_lim(1)*ysp+elec.chanpos(e,2), sz_lim(1)*zsp+elec.chanpos(e,3));
            % Color sphere
            % SHEILA: start from something like this: [~,cmap_ix] = min(abs(plot_dat(e,t_ix)-elec_cmap));
            set(e_sphr_ns{e}, 'EdgeColor', appropriate_color);%, 'FaceColor', ns_color, 'EdgeAlpha', edgealpha, 'FaceAlpha', facealpha);
        end
    end
    
    %SHEILA: update view_angle, potentially lighting too
    %view(view_angle(t_ix,:)); lighting gouraud;

    
    % Update time
    time_str.String = [num2str(hfa.time(t_ix)) ' s'];
    % Plot events
    if any([evnt_times{:}]==t_ix)
        for evnt_ix = 1:numel(plt_vars.evnt_types)
            if any(evnt_times{evnt_ix}==t_ix)
                evnt_str.String = plt_vars.evnt_types{evnt_ix};
                set(evnt_str,'Visible','on');
            end
        end
    else
        set(evnt_str,'Visible','off');
    end
    
    %pause(0.04);
    frames(t_ix) = getframe;
    
    % Remove all spheres
    %for all sphere; delete(e_sphr{e}); end;
end

%% Save movie
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/actv/'];
vid_out = VideoWriter([fig_dir plot_name vid_ext],vid_encoding);
vid_out.FrameRate = sample_rate*plt_vars.play_speed;
open(vid_out);
writeVideo(vid_out,frames);
close(vid_out);

end
