function fn_view_recon_stat_movie(SBJ, proc_id, stat_cond, atlas_id, roi_id, hemi, plot_out, plt_id, fig_vis, varargin)
%% Plot a reconstruction with electrodes colored according to statistics
%   PLT_VARS:
%       color = +/- or ROI
%       plot_nsig?
%       3 items below
%   FUTURE 2:
%       add option for stat_var to be a cell with 2nd stat for edge
%       gm_thresh [float] - threshold of GM % to include electrode (likely = 0)
%       z_thresh [float] - threshold of HFA z score to include electrode
%     future:
%       'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'PC': ANOVA of proportion congruence (red for sig)
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_cond [cell array] - {'stat_id','an_id','cond'}
%     stat_id [str] - ID of statistical analysis
%     an_id   [str] - analysis ID for preprocessing, filtering, etc.
%     cond    [str] - ID of condition/group to overlap with
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pCNI': ANOVA of previous trial congruence (red for sig)
%       'PC': ANOVA of proportion congruence (red for sig)
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%   hemi [str] - {'r','l','b'} for right, left, or both hemispheres (shouldn't be both here, too hard to see...)
%   plot_out [0/1] - plot electrodes outside of the hemisphere or ROIs of interest?
%   plt_id [str] - which set of plotting params
%   varargin [cell]:
%       view_angle [int,int] - angle to view the brain from
%       mesh_alpha [float] - transparency of the surface mesh

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

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

%% Prep variables
% Set up vars
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

cmap = eval([plt_vars.cmap_name '(64);']);
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

% Organize IDs
for st_ix = 1:numel(stat_cond)
    if numel(stat_cond)~=3; error('stat_cond must contain {stat_id, an_id, cond}'); end
end
stat_id = stat_cond{1};
an_id   = stat_cond{2};
cond_id = stat_cond{3};

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end

% Implement the default options
if ~exist('view_angle','var')
    if strcmp(hemi,'r')
        view_angle = [120 30];
    elseif strcmp(hemi,'l')
        view_angle = [-120 30];
    else
        error(['unknown hemi: ' hemi]);
    end
elseif strcmp(hemi,'b')
    error('no sense in making a both hemi movie!');
end
if ~exist('mesh_alpha','var')
    if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
        mesh_alpha = 0.4;
    else
        mesh_alpha = 0.8;
    end
end

%% Load elec struct
load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat']);

% Remove electrodes that aren't in hemisphere
if ~plot_out
    cfgs = [];
    cfgs.channel = fn_select_elec_lab_match(elec, hemi, atlas_id, []);
    elec = fn_select_elec(cfgs, elec);
end

%% Load brain recon
mesh = fn_load_recon_mesh(SBJ,'pat','','pial',hemi);

%% Load Stats
% Determine options: {'actv','CI','RT','CNI','PC'}
if strcmp(cond_id,'actv')
    % Load HFA and stat for colors and significance masking
    load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '.mat']);
    load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_' stat_id '.mat']);
    stat = actv;
    
    % Trim to channels and epoch
    cfgs = []; cfgs.latency = plt_vars.movie_lim;
    hfa  = ft_selectdata(cfgs,hfa);
%     elec = fn_select_elec(cfgs,elec);
    
    % Two sided color and size maps
    clim = [-max(abs(stat.avg(:))) max(abs(stat.avg(:)))];
    sz_map = [linspace(plt_vars.sz_lim(2), plt_vars.sz_lim(1), size(cmap,1)/2)...
              linspace(plt_vars.sz_lim(1), plt_vars.sz_lim(2), size(cmap,1)/2)];
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
%             elec_colors{ch_ix,1} = plt_vars.ns_color;
%         end
%     end    
elseif any(strcmp(cond_id,{'CNI','PC','pCNI','pCI'}))
    stat_fname = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id '_' an_id '.mat'];
    load(stat_fname,'w2','st');
    
    % Convert % explained variance to 0-100 scale
    w2.trial = w2.trial*100;
    stat = w2;
    cond_ix = strcmp(st.groups,cond_id);
    
    % Convert to ms timing
    win_times = 
    stat.mask = squeeze(stat.qval(cond_ix,:,:)<=st.alpha);
    
    % Get ROI colors
    elec.color = fn_roi2color(elec.(roi_field));
    
    % One sided color and size maps
    all_w2 = stat.trial(cond_ix,:,:);
%     clim = [prctile(all_w2(:),5) prctile(all_w2(:),95)];
    clim = [min(all_w2(:)) max(all_w2(:))];
    sz_map = linspace(plt_vars.sz_lim(1), plt_vars.sz_lim(2), size(cmap,1));
else
    error(['unknown stat_id: ' stat_id]);
end

%% Prepare stat plotting
% Confirm stat covers movie_lim
if any(plt_vars.movie_lim==-1)
    if plt_vars.movie_lim(1)==-1
        plt_vars.movie_lim(1) = stat.time(1);
    end
    if plt_vars.movie_lim(2)==-1
        plt_vars.movie_lim(2) = stat.time(end);
    end
end
if min(stat.time)>plt_vars.movie_lim(1)+0.0001 || max(stat.time)<plt_vars.movie_lim(2)-0.0001
    error('stat.time does not cover move_lim!');
end

% Get plotting time and colors
plot_time_idx = all([stat.time>=plt_vars.movie_lim(1); stat.time<=plt_vars.movie_lim(2)],1);
st_map = linspace(clim(1), clim(2), size(cmap,1));

%% Load timing info
evnt_time_idx = false([numel(plt_vars.evnt_type) numel(plot_time_idx)]);
for evnt_ix = 1:numel(plt_vars.evnt_type)
    if strcmp(plt_vars.evnt_type{evnt_ix},'S')
        [~,start_ix] = min(abs(stat.time-0));
        [~,end_ix]   = min(abs(stat.time-plt_vars.evnt_len));
    elseif strcmp(plt_vars.evnt_type{evnt_ix},'R')
        [~,start_ix] = min(abs(stat.time-mean(trial_info.response_time)));
        [~,end_ix]   = min(abs(stat.time-(mean(trial_info.response_time)+plt_vars.evnt_len)));
    end
    evnt_time_idx(evnt_ix,start_ix:end_ix) = true;
end
% Check events don't overlap
if size(evnt_time_idx,1)>1 && any(all(evnt_time_idx,1)); error('event label overlap!'); end

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
plot_name = [SBJ '_' cond_id '_' stat_id '_' an_id '_' hemi];
f = figure('Name',plot_name,'Visible',fig_vis); hold on;
frames = struct('cdata',[],'colormap',[]); frame_ix = 0;

% Plot 3D mesh
mesh_obj = ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
%SHEILA: this would be the first of potentially more than one view_angle
view(view_angle); material dull; lighting gouraud;
l = camlight;

[xsp, ysp, zsp] = sphere(250);

% Plot Timing info
% Text for time in sec (X,Y,Z in mm)
time_str = text(plt_vars.time_str_pos(1), plt_vars.time_str_pos(2), plt_vars.time_str_pos(3),...
    [num2str(stat.time(find(plot_time_idx,1)),'%.3f') ' s'], 'color', plt_vars.time_str_color,...
    'fontsize', plt_vars.time_str_size, 'horizontalalignment', plt_vars.time_str_horzalign);
evnt_str = text(plt_vars.evnt_str_pos(1), plt_vars.evnt_str_pos(2), plt_vars.evnt_str_pos(3),...
    plt_vars.evnt_type{1}, 'color', plt_vars.evnt_str_color,...
    'fontsize', plt_vars.evnt_str_size, 'horizontalalignment', plt_vars.evnt_str_horzalign,...
    'fontweight', plt_vars.evnt_str_weight);
set(evnt_str,'Visible','off');

if strcmp(cond_id,'actv')
    % Save base figure with colorbar
    cbar = colorbar; caxis(clim);
    fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/movies/' cond_id '/'];
    if ~exist(fig_dir,'dir')
        [~] = mkdir(fig_dir);
    end
    fig_fname = [fig_dir plot_name '.svg'];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
elseif any(strcmp(cond_id,{'CNI','PC','pCNI','pCI'})) && strcmp(plt_vars.evnt_type{1},'S')
    frame_ix = frame_ix + 1;
    % Plot electrodes at t = 0 (no actual stat for this, just for movie)
    e_sphr = gobjects(size(elec.label));
    for e = 1:numel(elec.label)
        ch_ix = find(strcmp(stat.label,elec.label{e}));
        % Create sphere
        sz = plt_vars.ns_sz;
        color = plt_vars.ns_color;
        % Draw and color electrodes
        e_sphr(e) = surf(sz*xsp+elec.chanpos(e,1), sz*ysp+elec.chanpos(e,2), sz*zsp+elec.chanpos(e,3));
        set(e_sphr(e), 'EdgeColor', color);%, 'FaceColor', plt_vars.ns_color, 'EdgeAlpha', edgealpha, 'FaceAlpha', facealpha);
    end
    
    % Grab t = 0 frame
    time_str.String = [num2str(0,'%.3f') ' s'];
    % Plot S event
    evnt_str.String = plt_vars.evnt_type{1};
    set(evnt_str,'Visible','on');
    
    %pause(plt_vars.frame_delay);
    frames(frame_ix) = getframe;
    
    % Remove all spheres to plot next round
    for e = 1:numel(elec.label)
        delete(e_sphr(e));
    end
end

%% Plot frames
fprintf('Starting frames (%i total):\n\t',sum(plot_time_idx)/plt_vars.frame_skip);
tic;
for t_ix = 1:plt_vars.frame_skip:numel(plot_time_idx)
    if plot_time_idx(t_ix)
        frame_ix = frame_ix + 1;
        if mod(frame_ix,5)==0; fprintf('%i..',t_ix); end
        if mod(frame_ix,100)==0; fprintf('\n\t'); end
        % For plotting times, plot electrodes
        e_sphr = gobjects(size(elec.label));
        for e = 1:numel(elec.label)
            ch_ix = find(strcmp(stat.label,elec.label{e}));
            % Create sphere
            if stat.mask(ch_ix,t_ix)
                if strcmp(cond_id,'actv')
                    [~,map_ix] = min(abs(st_map-stat.avg(ch_ix,t_ix)));
                    color = cmap(map_ix,:);
                else
                    [~,map_ix] = min(abs(st_map-stat.trial(cond_ix,ch_ix,t_ix)));
                    color = elec.color(e,:);
                end
                sz    = sz_map(map_ix);
            else
                sz = plt_vars.ns_sz;
                color = plt_vars.ns_color;
            end
            % Draw and color electrodes
            e_sphr(e) = surf(sz*xsp+elec.chanpos(e,1), sz*ysp+elec.chanpos(e,2), sz*zsp+elec.chanpos(e,3));
            set(e_sphr(e), 'EdgeColor', color);%, 'FaceColor', plt_vars.ns_color, 'EdgeAlpha', edgealpha, 'FaceAlpha', facealpha);
        end
        
        % Update view_angle and lighting
        %view(view_angle(t_ix,:)); lighting gouraud;
        
        % Update time
        time_str.String = [num2str(stat.time(t_ix),'%.3f') ' s'];
        % Plot events
        if any(evnt_time_idx(:,t_ix))
            evnt_str.String = plt_vars.evnt_type{evnt_time_idx(:,t_ix)~=0};
            set(evnt_str,'Visible','on');
        else
            set(evnt_str,'Visible','off');
        end
        
        %pause(plt_vars.frame_delay);
        frames(frame_ix) = getframe;
        
        % Remove all spheres to plot next round
        if t_ix~=numel(plot_time_idx)
            for e = 1:numel(elec.label)
                delete(e_sphr(e));
            end
        end
    end
end
fprintf('\nDONE: Elapsed time = %f\n',toc);

%% Save movie
frame_rate = 1/stat.time(2)-stat.time(1);
vid_out = VideoWriter([fig_dir plot_name plt_vars.vid_ext],plt_vars.vid_encoding);
vid_out.FrameRate = frame_rate*plt_vars.play_speed;
open(vid_out);
writeVideo(vid_out,frames);
close(vid_out);

end
