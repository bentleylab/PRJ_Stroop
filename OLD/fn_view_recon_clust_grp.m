function fn_view_recon_clust_grp(SBJs, clust_id, stat_id, an_id, atlas_id, roi_id,...
                    reg_type, show_labels, hemi, plot_out, plot_ns, plot_clusters, plot_recon,...
                    fig_vis, save_fig)%, view_angle)
                error('didnt finsih wrote reconatlas grp stt sphr clust instead')
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJ [str] - subject ID to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   plot_type [str] - {'ortho', '3d'} choose 3 slice orthogonal plot or 3D surface rendering
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   clust_id [str] - cluster analysis to pull
%   plot_out [0/1] - exclude electrodes that don't match atlas or aren't in hemisphere
%   plot_ns [0/1] - plot electrodes without significant effects?
%   plot_clusters [0/1] - plot time series (ANOVA and centroid) per bin
fig_filetype = 'svg';
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
view_space = 'mni';
pipeline_id = 'main_ft';

view_angle = [-90 0];
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end
% if strcmp(roi_id,'tissue') || strcmp(roi_id,'tissueC')
%     tis_suffix = '_tis';
% else
tis_suffix = '';
% end

%% Processing params
eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
if strcmp(stat_id,'actv') || strcmp(stat_id,'CSE')
    cond_lab = stat_id;
elseif strcmp(stat_id,'corrRT_CNI_pcon_WL200_WS50')
    % Get condition info
    [grp_lab, ~, ~] = fn_group_label_styles(model_lab);
    % if rt_correlation
    [rt_lab, ~, ~]     = fn_group_label_styles('RT');
    % end
    cond_lab = [grp_lab rt_lab];
else
    error(['Unknown stat_id: ' stat_id]);
end

if strcmp(an_id(1:5),'HGm_S')
    event_lab = 'stim';
    peak_bins = [0.3 0.6 1 2];
elseif strcmp(an_id(1:5),'HGm_R')
    event_lab = 'resp';
    peak_bins = [7 14 21 40; 7 14 21 40; 400 800 1200 2000];
end
stat_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% ROI info
[roi_list, ~] = fn_roi_label_styles(roi_id);
if strcmp(atlas_id,'Yeo7') || strcmp(atlas_id,'Yeo17')
    elec_space = 'mni_v';
else
    elec_space = 'pat';
end

% Cluster info
clust_names  = cell([1 numel(roi_list)]);
for clust_ix = 1:numel(roi_list)
    clust_names{clust_ix} = ['C' num2str(clust_ix)];
end
bin_colors = [0 0 1; 0 1 1; 1 0 1; 1 0 0];  %garish blue, cyan, magenta, red
% bin_colors = [27 158 119; 117 112 179; 217 95 2; 231 41 138]./255; %qualitative max diff from gROI, cool to hot
% bin_colors = [43 131 186; 171 221 164; 253 174 97; 215 25 28]./255;%diverging, green to red light middle
% bin_colors = [230 97 1; 253 184 99; 178 171 210; 94 60 153]./255;%diverging red to purple, too soft
% bin_colors = [240 249 232; 186 228 188; 123 204 196; 43 140 190]./255;
% bin_colors = [161, 218, 180; 65, 182, 196; 44 127 184; 37 52 148]./255;%darker, not enough contrast

%% Load cluster data
cond_ix = 1;
elec_sbj = cell([numel(SBJs) numel(cond_lab)]);
good_sbj = true([numel(SBJs) numel(cond_lab)]);
all_roi_labels = cell([numel(cond_lab) 1]);
all_roi_colors = cell([numel(cond_lab) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load elec struct (already in the clust .mat)
    % contains: 'clust_data','clusters','centroids','dist_sums','distances','elec'
    %   and maybe 'eva' if nCH criterion
    clust_fname = [SBJ_vars.dirs.proc SBJ '_' clust_id '_' stat_id '_' an_id '_' atlas_id '_' roi_id '.mat'];
    load(clust_fname);
    
    %% Load stats
    f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
    load(f_name);
    
    % Select analysis (ROI) elecs
    cfgs = []; cfgs.channel = elec.label;
    stat = ft_selectdata(cfgs,stat);
    w2 = ft_selectdata(cfgs,w2);
    
    % Get Sliding Window Parameters
    win_lim    = fn_sliding_window_lim(stat.time,win_len,win_step);
    win_center = round(mean(win_lim,2));
    
    % Convert % explained variance to 0-100 scale
    w2.trial = w2.trial*100;
    
    % Trim data to plotting epoch
    %   NOTE: stat should be on stat_lim(1):stat_lim(2)+0.001 time axis
    %   w2 should fit within that since it's averaging into a smaller window
    cfg_trim = [];
    if strcmp(event_lab,'stim')
        cfg_trim.latency = plt_vars.plt_lim_S;
    else
        cfg_trim.latency = plt_vars.plt_lim_R;
    end
    % hfa{1}  = ft_selectdata(cfg_trim,hfa{1});
    stat = ft_selectdata(cfg_trim,stat);
    sample_rate = (numel(stat.time)-1)/(stat.time(end)-stat.time(1));
    
    qvals = NaN(size(w2.pval));
    for ch_ix = 1:numel(stat.label)
        % FDR correct pvalues for ANOVA
        [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
    end
    
    %% Exclude non-sig elecs
    if ~plot_ns
        plot_ch = false([numel(elec.label) numel(cond_lab)]);
        for cond_ix = 1:numel(cond_lab)
            for ch_ix = 1:numel(elec.label)
                if any(strcmp(cond_lab{cond_ix},{'CNI','pcon'})) && any(squeeze(qvals(cond_ix,ch_ix,:))<0.05)
                    plot_ch(ch_ix,cond_ix) = true;
                elseif strcmp(cond_lab{cond_ix},'RT') && sum(squeeze(stat.mask(ch_ix,1,:)))>0
                    plot_ch(ch_ix,cond_ix) = true;
                end
            end
        end
    else
        plot_ch = true([numel(elec.label) numel(cond_lab)]);
    end
    
    %% Sort clusters by peak time
    clust_bin = cell([numel(cond_lab) 1]);
    for cond_ix = 1:numel(cond_lab)
        [~,peak_time] = max(centroids{cond_ix},[],2);
        [~,clust_bin{cond_ix}] = histc(peak_time,peak_bins(cond_ix,:));
    end
end

%% Plot centroids across SBJ
% Find plot limits
max_w2 = max(max(squeeze(w2.trial(1,:,:))));
min_w2 = min(min(squeeze(w2.trial(1,:,:))));
% max_w2 = max(max(max(w2.trial)));
% min_w2 = min(min(min(w2.trial)));
ylim1_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
ylims1  = [min_w2-ylim1_fudge max_w2+ylim1_fudge];
yticks1 = 0:1:ylims1(2);

max_rho = max(max(squeeze(stat.rho)));
min_rho = min(min(squeeze(stat.rho)));
ylims2  = [round(min_rho*10)/10-0.1 round(max_rho*10)/10+0.1]; % extra on top and bottom for StdErr
yticks2 = ylims2(1):0.1:ylims2(2);

if plot_clusters
    for cond_ix = 1%:numel(cond_lab)
        % Plot all ANOVA ts
        fig_name = [SBJ '_' cond_lab{cond_ix} '_' event_lab '_' roi_id '_' atlas_id '_binTS'];
        if ~plot_ns
            fig_name = [fig_name '_sig'];
        end
        f = figure('Name',fig_name,'Visible',fig_vis);
        ax = [];
        for ch_ix = 1:size(clusters,1)
            clust_n = clusters(ch_ix,cond_ix);
            time_bin = clust_bin{cond_ix}(clust_n)+1;
            subplot(4,1,time_bin); hold on;
            if plot_ch(ch_ix,cond_ix)
                if cond_ix <= numel(grp_lab)
                    plot(win_center,squeeze(w2.trial(cond_ix,ch_ix,:))',...
                        'Color',elec.roi_color(ch_ix,:));
                    % Find significant periods
                    sig_chunks = fn_find_chunks(squeeze(qvals(cond_ix,ch_ix,:))<0.05);
                    sig_chunks(squeeze(qvals(cond_ix,ch_ix,sig_chunks(:,1)))>0.05,:) = [];
                    for sig_ix = 1:size(sig_chunks,1)
                        line(win_center(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                            squeeze(w2.trial(cond_ix,ch_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                            'Color',elec.roi_color(ch_ix,:),...
                            'LineWidth',plt_vars.sig_width);
                    end
                else
                    % RT correlation significant time periods
                    plot(1:numel(stat.time),squeeze(stat.rho(ch_ix,:,:))',...
                        'Color',elec.roi_color(ch_ix,:),'LineStyle','-');
                    
                    sig_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,:,:)));
                    sig_chunks(squeeze(stat.mask(ch_ix,:,sig_chunks(:,1)))==0,:) = [];
                    for sig_ix = 1:size(sig_chunks,1)
                        sig_times = sig_chunks(sig_ix,1):sig_chunks(sig_ix,2);
                        line(sig_times,squeeze(stat.rho(ch_ix,:,sig_times)),...
                            'Color',elec.roi_color(ch_ix,:),'LineStyle',plt_vars.sig_style,...
                            'LineWidth',plt_vars.sig_width);
                    end
                end
            end
        end
        for clust_ix = 1:numel(clust_bin{cond_ix})            
            time_bin = clust_bin{cond_ix}(clust_ix)+1;
            subplot(4,1,time_bin);
            % Plot centroids for overall trends
            y_scale_factor = 0.15;
            if cond_ix <= numel(grp_lab)
                ylims = ylims1;
                yticks = yticks1;
                ylab = '% Variance Explained';
                y_data = centroids{cond_ix}(clust_ix,:);
                ylim_scaled = [ylims(1)+diff(ylims)*y_scale_factor ylims(2)-diff(ylims)*y_scale_factor];
                y_data = ylim_scaled(1) + [(y_data-min(y_data))./(max(y_data)-min(y_data))].*(ylim_scaled(2)-ylim_scaled(1));
                plot(win_center,y_data,...
                    'Color',bin_colors(time_bin,:),'LineWidth',2.5,'LineStyle',':');
            else
                ylims = ylims2;
                yticks = yticks2;
                ylab = 'Correlation with RT';
                y_data = centroids{cond_ix}(clust_ix,:);
                ylim_scaled = [ylims(1)+diff(ylims)*y_scale_factor ylims(2)-diff(ylims)*y_scale_factor];
                y_data = ylim_scaled(1) + [(y_data-min(y_data))./(max(y_data)-min(y_data))].*(ylim_scaled(2)-ylim_scaled(1));
                plot(1:numel(stat.time),y_data,...
                    'Color',bin_colors(time_bin,:),'LineWidth',2.5,'LineStyle',':');
            end
            
            % Plot event
            if strcmp(event_lab,'stim')
                x_tick_lab       = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
                mean_RT_pre = find(stat.time<=mean_RT);
                event_line = line([mean_RT_pre(end) mean_RT_pre(end)],ylims,...
                    'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                    'LineStyle',plt_vars.evnt_style);
            else
                x_tick_lab       = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
                event_line = line([find(stat.time==0) find(stat.time==0)],ylims,...
                    'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                    'LineStyle',plt_vars.evnt_style);
            end
            
            % Plotting parameters
%             ax(bin_ix).Box           = 'off';
            set(gca,'ylim',ylims);
            set(gca,'YTick',yticks);
%             set(gca,'YLabel',ylab);
            set(gca,'XLim',[0,size(stat.time,2)]);
            set(gca,'XTick',0:plt_vars.x_step_sz*sample_rate:size(stat.time,2));
            set(gca,'XTickLabel',x_tick_lab);
%             set(gca,'XLabel','Time (s)');
        end
        
        if save_fig
            fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/clust/' clust_id '/' stat_id '/' an_id '/'];
            fig_fname = [fig_dir fig_name '.' fig_filetype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
        end
%         % Plot centroids for overall trends
%         f = figure('Name',[cond_lab{cond_ix} '_centroid_' event_lab]);
%         ax = [];
%         for clust_ix = 1:size(centroids{cond_ix},1)
%             time_bin = clust_bin{cond_ix}(clust_ix)+1;
%             ax(bin_ix) = subplot(4,1,time_bin); hold on;
%             if cond_ix <= numel(grp_lab)
%                 plot(win_center,centroids{cond_ix}(clust_ix,:),'Color',bin_colors(time_bin,:));
%             else
%                 plot(1:numel(stat.time),centroids{cond_ix}(clust_ix,:),'Color',bin_colors(time_bin,:));
%             end
%             
%             % Plot event
%             if strcmp(event_lab,'stim')
%                 x_tick_lab       = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
%                 mean_RT_pre = find(stat.time<=mean_RT);
%                 event_line = line([mean_RT_pre(end) mean_RT_pre(end)],ylim,...
%                     'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
%                     'LineStyle',plt_vars.evnt_style);
%             else
%                 x_tick_lab       = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
%                 event_line = line([find(stat.time==0) find(stat.time==0)],ylim,...
%                     'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
%                     'LineStyle',plt_vars.evnt_style);
%             end
%             
%             % Plotting parameters
% %             ax(bin_ix).Box           = 'off';
% %             set(gca,'YLim',ylims);
% %             set(gca,'YTick',yticks);
% %             set(gca,'YLabel',ylab);
%             set(gca,'XLim',[0,size(stat.time,2)]);
%             set(gca,'XTick',0:plt_vars.x_step_sz*sample_rate:size(stat.time,2));
%             set(gca,'XTickLabel',x_tick_lab);
% %             set(gca,'XLabel','Time (s)');
%         end
    end
end

if plot_recon
    %% Remove electrodes that aren't in atlas ROIs
    if ~plot_out
        atlas_out_elecs = elec.label(strcmp(elec.atlas_label,'no_label_found'));
        if ~isempty(atlas_out_elecs)
            error('there shouldnt be elecs that arent in the ROIs, only hemi! check yoself...');
        end
        if ~strcmp(hemi,'b')
            plot_ch(~strcmp(elec.hemi,hemi),:) = false;
        end
    end
    
    %% Load brain recon
    if strcmp(view_space,'pat')
        if strcmp(hemi,'r') || strcmp(hemi,'l')
            mesh = ft_read_headshape([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_' hemi 'h.mat']);
        elseif strcmp(hemi,'b')
            mesh = ft_read_headshape({[SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_rh.mat'],...
                [SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_lh.mat']});
        else
            error(['Unknown hemisphere selected: ' hemi]);
        end
        mesh.coordsys = 'acpc';
    elseif strcmp(view_space,'mni')
        if strcmp(reg_type,'v')
            if strcmp(hemi,'r')
                load([ft_dir 'template/anatomy/surface_pial_right.mat']);
            elseif strcmp(hemi,'l')
                load([ft_dir 'template/anatomy/surface_pial_left.mat']);
            elseif strcmp(hemi,'b')
                load([ft_dir 'template/anatomy/surface_pial_both.mat']);
            else
                error(['Unknown hemisphere option: ' hemi]);
            end
            %         mesh.coordsys = 'mni';
        elseif strcmp(reg_type,'s')
            if strcmp(hemi,'r') || strcmp(hemi,'l')
                mesh = ft_read_headshape([root_dir 'PRJ_Stroop/data/atlases/freesurfer/fsaverage/' hemi 'h.pial']);
            elseif strcmp(hemi,'b')
                error('hemisphere "b" not yet implemented for reg_type: "srf"!');
                mesh = ft_read_headshape([ft_dir 'subjects/fsaverage/surf/' hemi 'h.pial']);
            else
                error(['Unknown hemisphere option: ' hemi]);
            end
            mesh.coordsys = 'fsaverage';
        else
            error(['Unknown registration type (reg_type): ' reg_type]);
        end
    else
        error(['Unknown view_space: ' view_space]);
    end
    
    %% Match elecs to atlas ROIs
    clust_n = unique(clusters);
    elec.clust = clusters;
    % clust_colors = distinguishable_colors(numel(clust_n));
    % if any(strcmp(atlas_id,{'DK','Dx','Yeo7'}))
    %     elec.roi       = fn_atlas2roi_labels(elec.atlas_label,atlas_id,roi_id);
    %     if strcmp(roi_id,'tissueC')
    %         elec.roi_color = fn_tissue2color(elec);
    %     elseif strcmp(atlas_id,'Yeo7')
    %         elec.roi_color = fn_atlas2color(atlas.name,elec.roi);
    %     else
    %         elec.roi_color = fn_roi2color(elec.roi);
    %     end
    % elseif any(strcmp(atlas_id,{'Yeo17'}))
    %     elec.roi       = elec.atlas_label;
    %     elec.roi_color = fn_atlas2color(atlas.name,elec.roi);
    % end
    
    %% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
    f = {};
    for cond_ix = 1%:numel(cond_lab)
        if strcmp(stat_id,'actv')
            plot_name = [SBJ '_actv_HFA_' an_id];
        elseif strcmp(stat_id,'CSE')
            plot_name = [SBJ '_' stat_id '_' an_id];
        else
            plot_name = [SBJ '_ANOVA_' cond_lab{cond_ix} '_' stat_id '_' an_id];
        end
        f{cond_ix} = figure('Name',plot_name,'Visible',fig_vis);
        
        % Plot 3D mesh
        mesh_alpha = 0.8;
        if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
            mesh_alpha = 0.3;
        end
        ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
        
        % Plot electrodes on top
        cfgs = [];
        for e = 1:numel(elec.label)
            if plot_ch(e,cond_ix)
                clust_n = clusters(e,cond_ix);
                time_bin = clust_bin{cond_ix}(clust_n)+1;
                cfgs.channel = elec.label{e};
                elec_tmp = fn_select_elec(cfgs, elec);
                if show_labels
                    ft_plot_sens(elec_tmp, 'elecshape', 'sphere', ...
                        'facecolor', squeeze(bin_colors(elec_tmp.clust(cond_ix),:)), 'label', 'label');
                else
                    ft_plot_sens(elec_tmp, 'elecshape', 'sphere', ...
                        'facecolor', squeeze(bin_colors(elec_tmp.clust(cond_ix),:)));
                end
            end
        end
        
        view(view_angle); material dull; lighting gouraud;
        l = camlight;
        fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
            'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
            '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
        set(f{cond_ix}, 'windowkeypressfcn',   @cb_keyboard);
    end
end

%% Plot SEEG data in 3D
% % Create volumetric mask of ROIs from fs parcellation/segmentation
% atlas = ft_read_atlas([SBJ_dir 'freesurfer/mri/aparc+aseg.mgz']);
% atlas.coordsys = 'acpc';
% cfg = [];
% cfg.inputcoord = 'acpc';
% cfg.atlas = atlas;
% cfg.roi = {'Right-Hippocampus', 'Right-Amygdala'};
% mask_rha = ft_volumelookup(cfg, atlas);

%% EXTRA CRAP:
% % Plot HFA from bipolar channels via clouds around electrode positions
% cfg = [];
% cfg.funparameter = 'powspctrm';
% cfg.funcolorlim = [-.5 .5];
% cfg.method = 'cloud';
% cfg.slice = '3d';
% cfg.nslices = 2;
% cfg.facealpha = .25;
% ft_sourceplot(cfg, freq_sel2, mesh_rha);
% view([120 40]); lighting gouraud; camlight;
%
% % 2D slice version:
% cfg.slice = '2d';
% ft_sourceplot(cfg, freq_sel2, mesh_rha);
%
% %% View grid activity on cortical mesh
% cfg = [];
% cfg.funparameter = 'powspctrm';
% cfg.funcolorlim = [-.5 .5];
% cfg.method = 'surface';
% cfg.interpmethod = 'sphere_weighteddistance';
% cfg.sphereradius = 8;
% cfg.camlight = 'no';
% ft_sourceplot(cfg, freq_sel, pial_lh);
%
% %% Prepare and plot 2D layout
% % Make layout
% cfg = [];
% cfg.headshape = pial_lh;
% cfg.projection = 'orthographic';
% cfg.channel = {'LPG*', 'LTG*'};
% cfg.viewpoint = 'left';
% cfg.mask = 'convex';
% cfg.boxchannel = {'LTG30', 'LTG31'};
% lay = ft_prepare_layout(cfg, freq);
% % Plot interactive
% cfg = [];
% cfg.layout = lay;
% cfg.showoutline = 'yes';
% ft_multiplotTFR(cfg, freq_blc);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)

if isempty(eventdata)
    % determine the key that corresponds to the uicontrol element that was activated
    key = get(h, 'userdata');
else
    % determine the key that was pressed on the keyboard
    key = parseKeyboardEvent(eventdata);
end
% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
    set(h, 'enable', 'off');
    drawnow;
    set(h, 'enable', 'on');
end

if strcmp(key, 'l') % reset the light position
    delete(findall(h,'Type','light')) % shut out the lights
    camlight; lighting gouraud; % add a new light from the current camera position
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = parseKeyboardEvent(eventdata)

key = eventdata.Key;
% handle possible numpad events (different for Windows and UNIX systems)
% NOTE: shift+numpad number does not work on UNIX, since the shift
% modifier is always sent for numpad events
if isunix()
    shiftInd = match_str(eventdata.Modifier, 'shift');
    if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
        % now we now it was a numpad keystroke (numeric character sent AND
        % shift modifier present)
        key = eventdata.Character;
        eventdata.Modifier(shiftInd) = []; % strip the shift modifier
    end
elseif ispc()
    if strfind(eventdata.Key, 'numpad')
        key = eventdata.Character;d
    end
end

if ~isempty(eventdata.Modifier)
  key = [eventdata.Modifier{1} '+' key];
end
