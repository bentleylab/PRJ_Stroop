function fn_view_recon_atlas_grp_ROI_stat_onset(SBJs, pipeline_id, stat_id, an_id, reg_type, show_labels,...
                                 hemi, atlas_id, roi_id, plot_roi, tbin_id, varargin)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJs [cell array str] - subject IDs to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of the stats
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pcon': ANOVA of proportion congruence (red for sig)
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the elecs
%       'gROI','mgROI','main3' - general ROIs (lobes or broad regions)
%       'ROI','thryROI','LPFC','MPFC','OFC','INS' - specific ROIs (within these larger regions)
%       'Yeo7','Yeo17' - colored by Yeo networks
%       'tissue','tissueC' - colored by tisseu compartment, e.g., GM vs WM vs OUT
%   plot_roi [str] - which surface mesh to plot

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Handle variables
% Error cases
if strcmp(hemi,'b') && ~strcmp(plot_roi,'OFC')
    error('hemi must be l or r for all non-OFC plots');
end
if ~any(strcmp(plot_roi,{'LPFC','MPFC','INS','OFC','TMP','PAR','lat','deep'}))
    error('plot_roi needs to be a lobe, "lat", or "deep"');
end

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'tick_step_s')
            tick_step_s = varargin{v+1};
        elseif strcmp(varargin{v},'view_angle')
            view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            mesh_alpha = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig')
            save_fig = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype')
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Define default options
if ~exist('tick_step_s','var')
    tick_step_s = 0.25;
end
if ~exist('view_angle','var')
    view_angle = fn_get_view_angle(hemi,plot_roi);
end
if ~exist('save_fig','var')
    save_fig = 0;
end
if ~exist('fig_ftype','var')
    fig_ftype = 'png';
end
% if ~exist('view_angle','var')
%     if strcmp(hemi,'l')
%         view_angle = [-60 -30];
%     elseif any(strcmp(hemi,{'r','b'}))
%         view_angle = [60 -30];
%     else
%         error(['unknown hemi: ' hemi]);
%     end
% end
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

%% Process stat_id
eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);
if strcmp(stat_id,'actv') || strcmp(stat_id,'CSE')
    cond_lab = stat_id;
    error('not ready for actv or CSE yet');
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

% Prep report
out_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP/' stat_id '/' an_id '/'];
sig_report_fname = [out_dir 'GRP_' atlas_id '_' plot_roi '_' event_type '_sig_report.txt'];
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');

%% Cluster info
if strcmp(event_type,'stim');
    error('onsets for stim need more thought!');
end
[roi_list, ~] = fn_roi_label_styles(roi_id);
if any(strcmp(plot_roi,{'deep','lat'}))
    [plot_roi_list, ~] = fn_roi_label_styles(plot_roi);
else
    plot_roi_list = {plot_roi};
end

% Get Time Bin and Sliding Window Parameters
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJs{1} '_vars.m']);
load(strcat(SBJ_vars.dirs.proc,SBJs{1},'_ANOVA_ROI_',stat_id,'_',an_id,'.mat'),'stat');
win_lim    = fn_sliding_window_lim(stat.time,win_len,win_step);
win_center = round(mean(win_lim,2));
% 4 ROIs = R time bins: -0.5, -0.1, 0.25, 0.6, 1.0
%   peak_bins = [7 14 21 40; 7 14 21 40; 400 800 1200 2000];
% 4 ROIs = R time bins: [0.3 0.6 1 2]
n_tbins    = cell(size(cond_lab));
bin_colors = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    % Determine coloring
    if strcmp(tbin_id,'cnts')
        if any(strcmp(cond_lab{cond_ix},grp_lab))
            n_tbins{cond_ix} = numel(win_center);
        else
            n_tbins{cond_ix} = numel(stat.time);
        end
    elseif strcmp(tbin_id,'eqROI')
        n_tbins{cond_ix} = numel(roi_list);
    else
        error('unknown tbin_id');
    end
    
    if n_tbins{cond_ix}==4 && ~strcmp(tbin_id,'cnts')
        bin_colors{cond_ix} = [0 0 1; 0 1 1; 1 0 1; 1 0 0];  %garish blue, cyan, magenta, red
    else
        bin_colors{cond_ix} = parula(n_tbins{cond_ix});
    end
end

%% Load Data
elec_sbj = cell([numel(SBJs) numel(cond_lab)]);
good_sbj = true([numel(SBJs) numel(cond_lab)]);
% Allocate fields that appendsens will remove
all_roi_labels = cell([numel(cond_lab) 1]);
all_roi_colors = cell([numel(cond_lab) 1]);
all_onset_ix   = cell([numel(cond_lab) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %% Elec
    % Load elec struct
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_mni',reg_suffix,'_',atlas_id,'_full.mat'];
    if exist([elec_fname(1:end-4) '_' roi_id '.mat'],'file')
        elec_fname = [elec_fname(1:end-4) '_' roi_id '.mat'];
    end
    tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
    
    % Append SBJ name to labels
    for elec_ix = 1:numel(elec_sbj{sbj_ix}.label)
        elec_sbj{sbj_ix}.label{elec_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix}.label{elec_ix}];
    end
    
    % Match elecs to atlas ROIs
    if any(strcmp(atlas_id,{'DK','Dx','Yeo7'}))
        if ~isfield(elec_sbj{sbj_ix,1},'man_adj')
            elec_sbj{sbj_ix,1}.roi       = fn_atlas2roi_labels(elec_sbj{sbj_ix,1}.atlas_lab,atlas_id,roi_id);
        end
        if strcmp(roi_id,'tissueC')
            elec_sbj{sbj_ix,1}.roi_color = fn_tissue2color(elec_sbj{sbj_ix,1});
        elseif strcmp(atlas_id,'Yeo7')
            elec_sbj{sbj_ix,1}.roi_color = fn_atlas2color(atlas_id,elec_sbj{sbj_ix,1}.roi);
        else
            elec_sbj{sbj_ix,1}.roi_color = fn_roi2color(elec_sbj{sbj_ix,1}.roi);
        end
    elseif any(strcmp(atlas_id,{'Yeo17'}))
        if ~isfield(elec_sbj{sbj_ix,1},'man_adj')
            elec_sbj{sbj_ix,1}.roi       = elec_sbj{sbj_ix,1}.atlas_lab;
        end
        elec_sbj{sbj_ix,1}.roi_color = fn_atlas2color(atlas_id,elec_sbj{sbj_ix,1}.roi);
    end
    
    % Select elecs matching hemi, atlas, and plot_roi_list
    plot_elecs = zeros([numel(elec_sbj{sbj_ix,1}.label) numel(plot_roi_list)]);
    for roi_ix = 1:numel(plot_roi_list)
        plot_elecs(:,roi_ix) = strcmp(elec_sbj{sbj_ix,1}.roi,plot_roi_list{roi_ix});
    end
    cfgs = [];
    cfgs.channel = intersect(fn_select_elec_lab_match(elec_sbj{sbj_ix,1}, hemi, atlas_id, roi_id),...
                             elec_sbj{sbj_ix,1}.label(any(plot_elecs,2)));
    elec_sbj{sbj_ix,1} = fn_select_elec(cfgs, elec_sbj{sbj_ix,1});
    
    % Create dummy onset variable
    elec_sbj{sbj_ix,1}.onset    = NaN(size(elec_sbj{sbj_ix,1}.label));
    elec_sbj{sbj_ix,1}.onset_ix = NaN(size(elec_sbj{sbj_ix,1}.label));
    % Copy for other conditions
    for cond_ix = 2:numel(cond_lab)
        elec_sbj{sbj_ix,cond_ix} = elec_sbj{sbj_ix,1};
    end
    
    %% Stats
    if isempty(elec_sbj{sbj_ix,1}.label)
        fprintf('%s has no elecs in %s, skipping...\n',SBJ,plot_roi);
        for cond_ix = 1:numel(cond_lab)
            elec_sbj{sbj_ix,cond_ix} = {};
            good_sbj(sbj_ix,cond_ix) = false;
        end
    else
        % Load data
        load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id,'.mat'));
        
        % FDR correct pvalues for ANOVA
        qvals = NaN(size(w2.pval));
        for elec_ix = 1:numel(stat.label)
            [~, ~, ~, qvals(:,elec_ix,:)] = fdr_bh(squeeze(w2.pval(:,elec_ix,:)));%,0.05,'pdep','yes');
        end
        
        % Aggregate results per elec
        for elec_ix = 1:numel(elec_sbj{sbj_ix,1}.label)
            orig_lab = strrep(elec_sbj{sbj_ix,1}.label{elec_ix},[SBJs{sbj_ix} '_'],'');
            stat_ix  = find(strcmp(stat.label,orig_lab));
            % Get ANOVA group onsets
            for grp_ix = 1:numel(grp_lab)
                if any(squeeze(qvals(grp_ix,stat_ix,:))<0.05)
                    %                 sig_ch{sbj_ix,grp_ix} = [sig_ch{sbj_ix,grp_ix} stat.label{ch_ix}];
                    sig_onset_ix = find(squeeze(qvals(grp_ix,stat_ix,:))<0.05,1);
                    sig_onsets   = stat.time(win_lim(sig_onset_ix,1));
                    if strcmp(event_type,'resp')
                        elec_sbj{sbj_ix,grp_ix}.onset_ix(elec_ix) = sig_onset_ix(1);
                        elec_sbj{sbj_ix,grp_ix}.onset(elec_ix)    = sig_onsets(1);
                    elseif strcmp(event_type,'stim') && (sig_onsets(1)<mean_RTs(sbj_ix))
                        elec_sbj{sbj_ix,grp_ix}.onset_ix(elec_ix) = sig_onset_ix(1);
                        elec_sbj{sbj_ix,grp_ix}.onset(elec_ix)    = sig_onsets(1);
                    end
                end
            end
            
            % Get RT correlation onset
            if sum(squeeze(stat.mask(stat_ix,1,:)))>0
                %             sig_ch{sbj_ix,numel(grp_lab)+1} = [sig_ch{sbj_ix,numel(grp_lab)+1} stat.label{ch_ix}];
                mask_chunks = fn_find_chunks(squeeze(stat.mask(stat_ix,1,:)));
                mask_chunks(squeeze(stat.mask(stat_ix,1,mask_chunks(:,1)))==0,:) = [];
                % Convert the first onset of significance to time
                onset_time = stat.time(mask_chunks(1,1));
                % Exclude differences after the mean RT for this SBJ
                if strcmp(event_type,'resp')
                    elec_sbj{sbj_ix,numel(grp_lab)+1}.onset_ix(elec_ix) = mask_chunks(1,1);
                    elec_sbj{sbj_ix,numel(grp_lab)+1}.onset(elec_ix)    = onset_time;
                elseif strcmp(event_type,'stim') && (onset_time<mean_RTs(sbj_ix))
                    elec_sbj{sbj_ix,numel(grp_lab)+1}.onset_ix(elec_ix) = mask_chunks(1,1);
                    elec_sbj{sbj_ix,numel(grp_lab)+1}.onset(elec_ix)    = onset_time;
                end
            end
        end
        
        %     % Normalize all onset times by mean reaction time
        %     if strcmp(event_type,'stim')
        %         for roi_ix = 1:size(all_onsets,2)
        %             for cond_ix = 1:size(all_onsets,3)
        %                 all_onsets{sbj_ix,roi_ix,cond_ix} = all_onsets{sbj_ix,roi_ix,cond_ix}./mean_RTs(sbj_ix);
        %             end
        %         end
        %     end
        
        %     % Determine options: {'actv','CI','RT','CNI','pcon'}
        %     sig_ch = cell(size(cond_lab));
        %     if strcmp(stat_id,'actv')
        %         load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_actv_mn100.mat'],'actv_ch');
        %         sig_ch{1} = actv_ch;
        %         clear actv_ch
        %     elseif strcmp(stat_id,'CSE')
        %         load([SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'_',stat_id,'.mat'],'stat');
        %         for elec_ix = 1:numel(stat.label)
        %             if any(stat.mask(elec_ix,1,:))
        %                 sig_ch{1} = [sig_ch{1} stat.label(elec_ix)];
        %             end
        %         end
        %         clear stat
        %     else    % ANOVA
        %         eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
        %         f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
        %         load(f_name,'stat','w2');
        %
        %         % FDR correct pvalues for ANOVA
        %         for elec_ix = 1:numel(stat.label)
        %             pvals = squeeze(w2.pval(:,elec_ix,:));
        %             [~, ~, ~, qvals] = fdr_bh(pvals);%,0.05,'pdep','yes');
        %
        %             % Consolidate to binary sig/non-sig
        %             for cond_ix = 1:numel(cond_lab)
        %                 if strcmp(cond_lab{cond_ix},'RT') && any(stat.mask(elec_ix,1,:))
        %                     sig_ch{cond_ix} = [sig_ch{cond_ix} {[SBJs{sbj_ix} '_' stat.label{elec_ix}]}];
        %                 elseif any(strcmp(cond_lab{cond_ix},{'CNI','pcon'})) && any(qvals(cond_ix,:)<0.05,2)
        %                     sig_ch{cond_ix} = [sig_ch{cond_ix} {[SBJs{sbj_ix} '_' w2.label{elec_ix}]}];
        %                 end
        %             end
        %         end
        %         clear stat w2
        %     end
        
        % Select sig elecs
        fprintf(sig_report,'===============================================================\n');
        for cond_ix = 1:numel(cond_lab)
            sig_ix = find(~isnan(elec_sbj{sbj_ix,cond_ix}.onset));
            % Report overall numbers on significant electrodes for this SBJ
            fprintf(sig_report,'----------\t%s\t----------\n',cond_lab{cond_ix});
            fprintf(sig_report,'\t%s - %s  = %i / %i (%.02f) sig elecs:\n',SBJs{sbj_ix},cond_lab{cond_ix},...
                numel(sig_ix),numel(elec_sbj{sbj_ix,cond_ix}.label),...
                100*numel(sig_ix)/numel(elec_sbj{sbj_ix,cond_ix}.label));
            % Report individual onsets
            for ch_ix = 1:numel(sig_ix)
                elec_ix = sig_ix(ch_ix);%find(strcmp(elec_sbj{sbj_ix,cond_ix}.label,sig_ix{cond_ix}{sig_ix}));
                fprintf(sig_report,'%s - %s (%s) = %.3f\n',elec_sbj{sbj_ix,cond_ix}.label{elec_ix},elec_sbj{sbj_ix,cond_ix}.atlas_lab{elec_ix},...
                    elec_sbj{sbj_ix,cond_ix}.hemi{elec_ix}, elec_sbj{sbj_ix,cond_ix}.onset(elec_ix));
            end
            
            % Select sig elecs && elecs matching atlas
            if isempty(sig_ix)   % fn_select_elec errors if no elec left
                elec_sbj{sbj_ix,cond_ix} = {};
                good_sbj(sbj_ix,cond_ix) = false;
                warning([SBJ ' ' cond_lab{cond_ix} ' EMPTY!']);
                fprintf(sig_report,'\t!!! 0 sig_ch remain\n');
            else
                if numel(elec_sbj{sbj_ix,cond_ix}.label)>1
                    % Remove any remaining non-sig elecs
                    cfgs = [];
                    cfgs.channel = elec_sbj{sbj_ix,cond_ix}.label(sig_ix);
                    elec_sbj{sbj_ix,cond_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix,cond_ix});
                end
                all_roi_labels{cond_ix} = [all_roi_labels{cond_ix}; elec_sbj{sbj_ix,cond_ix}.roi];
                all_roi_colors{cond_ix} = [all_roi_colors{cond_ix}; elec_sbj{sbj_ix,cond_ix}.roi_color];
                all_onset_ix{cond_ix}     = [all_onset_ix{cond_ix}; elec_sbj{sbj_ix,cond_ix}.onset_ix];
                fprintf(sig_report,'\t%i sig_ch remain\n',numel(elec_sbj{sbj_ix,cond_ix}.label));
            end
            fprintf(sig_report,'---------------------------------------------------------------\n');
        end
    end
    fprintf(sig_report,'===============================================================\n');
    clear SBJ SBJ_vars SBJ_vars_cmd
end

%% Combine elec structs
elec = cell([numel(cond_lab) 1]);
for cond_ix = 1:numel(cond_lab)
    if any(good_sbj(:,cond_ix))
        elec{cond_ix} = ft_appendsens([],elec_sbj{good_sbj(:,cond_ix),cond_ix});
        % appendsens strips these fields
        elec{cond_ix}.roi       = all_roi_labels{cond_ix};
        elec{cond_ix}.roi_color = all_roi_colors{cond_ix};
        elec{cond_ix}.onset_ix  = all_onset_ix{cond_ix};
    else
        fprintf(2,'WARNING!!! No SBJ had sig elecs in %s for %s!\n',plot_roi,cond_lab{cond_ix});
    end
end

%% Load Atlas
atlas = fn_load_recon_atlas([],atlas_id);

% Get Atlas-ROI mapping
atlas_labels = fn_atlas_roi_select_mesh(atlas_id, plot_roi, hemi);
if strcmp(plot_roi,'deep')
    mtl_ix = ~cellfun(@isempty,strfind(atlas_labels,'Hippocampus')) | ...
             ~cellfun(@isempty,strfind(atlas_labels,'Amygdala'));
    mtl_labels = atlas_labels(mtl_ix);
    atlas_labels = atlas_labels(~mtl_ix);
end

%% Select ROI mesh
cfg = [];
cfg.inputcoord = atlas.coordsys;
cfg.atlas = atlas;
cfg.roi = atlas_labels;
roi_mask = ft_volumelookup(cfg,atlas);
seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = roi_mask;

if exist('mtl_labels','var')
    cfg.roi = mtl_labels;
    mtl_mask = ft_volumelookup(cfg,atlas);
    
    mtl_seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    mtl_seg.brain = mtl_mask;
end

cfg = [];
cfg.method      = 'iso2mesh';   % surface toolbox Arjen found
cfg.radbound    = 2;            % scalar indicating the radius of the target surface mesh element bounding sphere
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 100000;
cfg.smooth      = 3;
roi_mesh = ft_prepare_mesh(cfg, seg);
if exist('mtl_seg','var')
    mtl_mesh = ft_prepare_mesh(cfg, mtl_seg);
end

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
f = cell(size(cond_lab));
for cond_ix = 1%:numel(cond_lab)
    if strcmp(stat_id,'actv') || strcmp(stat_id,'CSE')
        fig_name = ['GRP_' stat_id '_' event_type '_' atlas_id '_' roi_id '_' plot_roi '_' hemi];
    else
        fig_name = ['GRP_' cond_lab{cond_ix} '_' event_type '_' atlas_id '_' roi_id '_' plot_roi '_' hemi];
    end
    f{cond_ix} = figure('Name',fig_name);
    
    % Plot 3D mesh
    ft_plot_mesh(roi_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
    if exist('mtl_mesh','var')
        ft_plot_mesh(mtl_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
    end

    % Plot electrodes on top
    if any(good_sbj(:,cond_ix))
        if numel(elec{cond_ix}.label)>1
            for e = 1:numel(elec{cond_ix}.label)
                cfgs = []; cfgs.channel = elec{cond_ix}.label(e);
                elec_tmp = fn_select_elec(cfgs,elec{cond_ix});
                ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', ...
                    bin_colors{cond_ix}(elec_tmp.onset_ix,:), 'label', lab_arg);
            end
        else
            ft_plot_sens(elec{cond_ix}, 'elecshape', 'sphere', 'facecolor', ...
                bin_colors{cond_ix}(elec{cond_ix}.onset_ix,:), 'label', lab_arg);
        end
    end
    
    % Colorbar
    caxis([min(stat.time) max(stat.time)]);
    colorbar('Ticks',min(stat.time):tick_step_s:max(stat.time));%'TickLabels',{'Cold','hot'}
    
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
    fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
        'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
        '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
    set(f{cond_ix}, 'windowkeypressfcn',   @cb_keyboard);
    
    if save_fig
        fig_fname = [out_dir fig_name fig_ftype];
        saveas(gcf,fig_fname);
    end
end

