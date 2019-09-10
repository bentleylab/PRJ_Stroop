function fn_view_recon_stat_movie(SBJ, proc_id, stat_cond, atlas_id, roi_id, hemi, mirror, plot_out, plt_id, fig_vis, varargin)
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
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pCNI': ANOVA of previous trial congruence (red for sig)
%       'PC': ANOVA of proportion congruence (red for sig)
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%   hemi [str] - {'r','l'} for right or left hemispheres ('b' gives error)
%   mirror [0/1] - plot the other hemi, 
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
proc_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Organize IDs
for st_ix = 1:numel(stat_cond)
    if numel(stat_cond)~=3; error('stat_cond must contain {stat_id, an_id, cond}'); end
end
stat_id = stat_cond{1};
an_id   = stat_cond{2};
cond_id = stat_cond{3};

%% Implement default input options
if ~exist('view_angle','var')
    if strcmp(hemi,'r')
        view_angle = [120 10];
    elseif strcmp(hemi,'l')
        view_angle = [-120 10];
    elseif strcmp(hemi,'b')
        error('no sense in making a both hemi movie! (use mirroring...)');
    else
        error(['unknown hemi: ' hemi]);
    end
    view_str = 'def';
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
    if mirror
        cfgs.channel = fn_select_elec_lab_match(elec, 'b', atlas_id, []);
    else
        cfgs.channel = fn_select_elec_lab_match(elec, hemi, atlas_id, []);
    end
    elec = fn_select_elec(cfgs, elec);
end

% Mirror hemispheres
if mirror
    elec.chanpos(~strcmp(elec.hemi,hemi),1) = -elec.chanpos(~strcmp(elec.hemi,hemi),1);
    hemi_str = [hemi 'b'];
end

%% Organize Movie Settings
% Create time index
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
% Adjust movie_lim to special cases
for lim_ix = 1:2
    if ~isnumeric(plt_vars.movie_lim{lim_ix})
        if strcmp(plt_vars.movie_lim{lim_ix},'D')
            plt_vars.movie_lim{lim_ix} = median(trial_info.response_time);
        end
    end
end
stat.time = plt_vars.movie_lim{1}:1/proc.resample_freq:plt_vars.movie_lim{2};

% Adjust view angle if custom
if ischar(view_angle)
    view_str = view_angle;
    if (strcmp(hemi,'l') && strcmp(view_angle,'med')) || (strcmp(hemi,'r') && strcmp(view_angle,'lat'))
        view_angle = [90 0];
    elseif (strcmp(hemi,'l') && strcmp(view_angle,'lat')) || (strcmp(hemi,'r') && strcmp(view_angle,'med'))
        view_angle = [-90 0];
    end
end

%% Load Stats
% Determine options: {'actv','CI','RT','CNI','PC'}
if strcmp(cond_id,'actv')
    % Load HFA and stat for colors and significance masking
%     load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '.mat']);
    load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_' stat_id '.mat']);
%     if numel(actv.label)~=numel(hfa.label) || ~all(strcmp(actv.label,hfa.label))
%         error('mismatched elecs in hfa and actv!');
%     end
    
    % Convert actv to stat struct
    stat.data = nan([numel(elec.label) numel(stat.time)]);
    warn = 0;
    for e_ix = 1:numel(elec.label)
        ch_ix = find(strcmp(actv.label,elec.label{e_ix}));
        for t_ix = 1:numel(stat.time)
            [~,actv_ix] = min(abs(actv.time-stat.time(t_ix)));
            if ~isnan(stat.data(e_ix,t_ix)) && ~warn; warning('time overlap!'); end
            if actv.mask(ch_ix,actv_ix)
                stat.data(e_ix,t_ix) = actv.avg(ch_ix,actv_ix);
            end
        end
    end
    
    % Define event timing
    if strcmp(st.evnt_lab,'R')
        st.ep_lab = 'DR';
        stat.evnt = ones(size(stat.time));
        for t_ix = 1:numel(stat.time)
            if stat.time(t_ix)>=0
                stat.evnt(t_ix) = 2;
            end
        end
    else
        st.ep_lab = 'BSDR'; % force inclusion of all possible events
        stat.evnt = zeros(size(stat.time)); warn = 0;
        for t_ix = 1:numel(stat.time)
            if stat.time(t_ix)<0
                % B: before t = 0
                stat.evnt(t_ix) = 1;
            elseif stat.time(t_ix)<min(trial_info.response_time)
                % S: 0 to min(RT)
                stat.evnt(t_ix) = 2;
            elseif stat.time(t_ix)<prctile(trial_info.response_time,75)
                % D: min(RT) to 75th percentile of RTs
                stat.evnt(t_ix) = 3;
            else
                % R: 75th percentile of RTs to end
                stat.evnt(t_ix) = 4;
            end
        end
    end
    
    % Two sided color and size maps
    clim = [-max(abs(actv.avg(:))) max(abs(actv.avg(:)))];
    sz_map = [linspace(plt_vars.sz_lim(2), plt_vars.sz_lim(1), size(plt_vars.cmap,1)/2)...
              linspace(plt_vars.sz_lim(1), plt_vars.sz_lim(2), size(plt_vars.cmap,1)/2)];
    
elseif any(strcmp(cond_id,{'CNI','PC','pCNI','pCI'}))
    stat_fname = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id '_' an_id '.mat'];
    load(stat_fname,'w2','st');
    cond_ix = find(strcmp(st.groups,cond_id));
    
    % Convert w2 to stat struct
    stat.data = nan([numel(elec.label) numel(stat.time)]);
    for e_ix = 1:numel(elec.label)
        ch_ix = find(strcmp(w2.label,elec.label{e_ix}));
        if any(w2.qval(cond_ix,ch_ix,:)<=st.alpha)
            sig_mat = nan([numel(w2.time) numel(stat.time)]);
            win_mat = nan([numel(w2.time) numel(stat.time)]);
            for w_ix = 1:numel(w2.time)
                % Grab data and mask (convert % explained variance to 0-100 scale)
                win_idx = find(stat.time>=w2.win_lim_s(w_ix,1) & stat.time<=w2.win_lim_s(w_ix,2));
                win_mat(w_ix,win_idx) = w2.trial(cond_ix,ch_ix,w_ix)*100;
                sig_mat(w_ix,win_idx) = w2.qval(cond_ix,ch_ix,w_ix)<=st.alpha;
            end
            % Average w2 across windows to single time series
            stat.data(e_ix,:) = nanmean(win_mat,1);
            % Remove data from epochs not overlapping with sig windows
            stat.data(e_ix,~any(sig_mat,1)) = nan;
        end
    end
    
    % Define event timing
    if isfield(w2,'stat_id_ts')
        stat.evnt = zeros(size(stat.time)); warn = 0;
        for w_ix = 1:numel(w2.time)
            win_idx = find(stat.time>=w2.win_lim_s(w_ix,1) & stat.time<=w2.win_lim_s(w_ix,2));
            if any(stat.evnt(win_idx)~=0) && ~warn; warning('evnt overlap!'); warn = 1; end
            stat.evnt(win_idx) = w2.stat_id_ts(w_ix);
        end
    else
        stat.evnt = ones(size(stat.time));
    end
    
    % Get ROI colors
    if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep'}))
        elec.color = fn_roi2color(elec.gROI);
    else
        elec.color = fn_roi2color(elec.ROI);
    end
    
    % One sided color and size maps
    all_w2 = w2.trial(cond_ix,:,:)*100;
%     clim = [prctile(all_w2(:),5) prctile(all_w2(:),95)];
    clim = [min(all_w2(:)) max(all_w2(:))];
    sz_map = linspace(plt_vars.sz_lim(1), plt_vars.sz_lim(2), size(plt_vars.cmap,1));
else
    error(['unknown stat_id: ' stat_id]);
end

% Get plotting time and colors
st_map = linspace(clim(1), clim(2), size(plt_vars.cmap,1));
% if min(stat.time)>plt_vars.movie_lim(1)+0.0001 || max(stat.time)<plt_vars.movie_lim(2)-0.0001
%     error('stat.time does not cover move_lim!');
% end

%% Load brain recon
mesh = fn_load_recon_mesh(SBJ,'pat','','pial',hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
plot_name = [SBJ '_' cond_id '_' stat_id '_' an_id '_' hemi_str '_' view_str];
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/movies/' cond_id '/'];
if ~exist(fig_dir,'dir')
    [~] = mkdir(fig_dir);
end

f = figure('Name',plot_name,'Visible',fig_vis); hold on;
frames = struct('cdata',[],'colormap',[]); frame_ix = 0;

% Plot 3D mesh
mesh_obj = ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
view(view_angle); material dull; lighting gouraud; camlight;

% Plot elecs
[xsp, ysp, zsp] = sphere(plt_vars.sphere_sz);
e_sphr = gobjects(size(elec.label));
for e = 1:numel(elec.label)
    % Draw and color electrodes
    sz = plt_vars.ns_sz;
    color = plt_vars.ns_color;
    e_sphr(e) = surf(sz*xsp+elec.chanpos(e,1), sz*ysp+elec.chanpos(e,2), sz*zsp+elec.chanpos(e,3));
    set(e_sphr(e), 'EdgeColor', color);%, 'FaceColor', plt_vars.ns_color, 'EdgeAlpha', edgealpha, 'FaceAlpha', facealpha);
end

% Plot Timing info
% Text for time in sec (X,Y,Z in mm)
time_str = text(plt_vars.time_str_pos(1), plt_vars.time_str_pos(2), plt_vars.time_str_pos(3),...
    [num2str(stat.time(1),'%.3f') ' s'], 'color', plt_vars.time_str_color,...
    'fontsize', plt_vars.time_str_size, 'horizontalalignment', plt_vars.time_str_horzalign);
evnt_str = text(plt_vars.evnt_str_pos(1), plt_vars.evnt_str_pos(2), plt_vars.evnt_str_pos(3),...
    st.ep_lab(stat.evnt(1)), 'color', plt_vars.evnt_str_color,...
    'fontsize', plt_vars.evnt_str_size, 'horizontalalignment', plt_vars.evnt_str_horzalign,...
    'fontweight', plt_vars.evnt_str_weight);

if strcmp(cond_id,'actv')
    % Save base figure with colorbar
    cbar = colorbar; caxis(clim);
    cbar_dir = [fig_dir 'cbar_snap/'];
    if ~exist(cbar_dir,'dir')
        [~] = mkdir(cbar_dir);
    end
    fig_fname = [cbar_dir plot_name '.svg'];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot frames
fprintf('Starting frames (%i total):\n',numel(stat.time)/plt_vars.frame_skip);
frame_plot_times = zeros(size(stat.time));
frame_plot_dur   = zeros(size(stat.time));
update_idx       = zeros(size(stat.time));
tic;
for t_ix = 1:plt_vars.frame_skip:numel(stat.time)
    frame_ix = frame_ix + 1;
%     if mod(frame_ix,5)==0; fprintf('%i..',t_ix); end
%     if mod(frame_ix,100)==0; fprintf('\n\t'); end
    % Draw new elecs if first frame or change in data
    if t_ix==1 || ~isequaln(stat.data(:,t_ix),stat.data(:,t_ix-1)) % ~= doesn't handle nans
        % Remove all spheres to plot next round
        % delete(e_sphr); e_sphr = gobjects(size(elec.label));
        for e = 1:numel(elec.label)
            if t_ix==1 || ~isequaln(stat.data(e,t_ix),stat.data(e,t_ix-1))
                update_idx(t_ix) = update_idx(t_ix) + 1;
                % Create sphere
                if ~isnan(stat.data(e,t_ix))
                    [~,map_ix] = min(abs(st_map-stat.data(e,t_ix)));
                    sz    = sz_map(map_ix);
                    if strcmp(cond_id,'actv')
                        color = plt_vars.cmap(map_ix,:);
                    else
                        color = elec.color(e,:);
                    end
                else
                    sz = plt_vars.ns_sz;
                    color = plt_vars.ns_color;
                end
                % Draw and color electrodes
                e_sphr(e) = surf(sz*xsp+elec.chanpos(e,1), sz*ysp+elec.chanpos(e,2), sz*zsp+elec.chanpos(e,3));
                set(e_sphr(e), 'EdgeColor', color, 'FaceColor', color);%, 'EdgeAlpha', edgealpha, 'FaceAlpha', facealpha);
            end
        end
    end
    % Update view_angle and lighting
    % view(view_angle(t_ix,:)); lighting gouraud;
    % delete(findall(h,'Type','light')) % shut out the lights
    % camlight; lighting gouraud; % add a new light from the current camera position

    % Update time and event strings
    time_str.String = [num2str(stat.time(t_ix),'%.3f') ' s'];
    evnt_str.String = st.ep_lab(stat.evnt(t_ix));
    
    %pause(plt_vars.frame_delay);
    frames(frame_ix) = getframe;
    frame_plot_times(t_ix) = toc;
    if t_ix==1
        frame_plot_dur(t_ix) = frame_plot_times(t_ix);
        fprintf('\t%i (%.1f, %i) = %.1f\n',t_ix,100*t_ix/numel(stat.time),...
            update_idx(t_ix),frame_plot_dur(t_ix));
    else
        frame_plot_dur(t_ix) = frame_plot_times(t_ix)-frame_plot_times(t_ix-1);
        fprintf('\t%i (%.1f, %i) = %.1f (%.3f)\n',t_ix,100*t_ix/numel(stat.time),...
            update_idx(t_ix),frame_plot_times(t_ix),frame_plot_dur(t_ix));
    end
end

%% Save movie
fprintf('Saving %s...\n',[fig_dir plot_name plt_vars.vid_ext]);
vid_out = VideoWriter([fig_dir plot_name plt_vars.vid_ext],plt_vars.vid_encoding);
vid_out.FrameRate = proc.resample_freq*plt_vars.play_speed;
open(vid_out);
writeVideo(vid_out,frames);
close(vid_out);

%% Print timing stats
fprintf('\nDONE: Elapsed time = %f\n',toc);
fprintf('frame_times (mean, non-update, update - SD) = (%.3f, %.3f, %.3f - %.3f)\n',...
    mean(frame_plot_dur), mean(frame_plot_dur(update_idx==0)),...
    mean(frame_plot_dur(update_idx>0)), std(frame_plot_dur(update_idx>0)));
fprintf('updates (n, min, mean - SD, max) = (%i, %i, %.2f - %.3f, %i)\n',sum(update_idx>0),min(update_idx(update_idx>0)),...
    mean(update_idx(update_idx>0)), std(update_idx(update_idx>0)), max(update_idx(update_idx>0)));
figure;scatter(frame_plot_dur,update_idx);
xlabel('time (s)'); ylabel('# elecs updated');

end
