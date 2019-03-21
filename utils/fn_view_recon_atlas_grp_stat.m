function fn_view_recon_atlas_grp_stat(SBJs, pipeline_id, stat_id, an_id, reg_type, show_labels,...
                                 hemi, atlas_id, roi_id, plot_out, varargin)
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

% Define default options
% view_space = 'mni';
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
fprintf('Using atlas: %s\n',atlas_id);

%% Process stat_id
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

% Prep report
out_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_reports/'];
sig_report_fname = [out_dir 'GRP_' stat_id '_' an_id '_' atlas_id '_' roi_id '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');

%% Load Data
elec_sbj = cell([numel(SBJs) numel(cond_lab)]);
good_sbj = true([numel(SBJs) numel(cond_lab)]);
all_roi_labels = cell([numel(cond_lab) 1]);
all_roi_colors = cell([numel(cond_lab) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load elec struct
    try
        elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_mni',reg_suffix,'_',atlas_id,'_full.mat'];
        if exist([elec_fname(1:end-4) '_' roi_id '.mat'],'file')
            elec_fname = [elec_fname(1:end-4) '_' roi_id '.mat'];
        end
        tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
    catch
        error([elec_fname 'doesnt exist, exiting...']);
    end
    
    % Append SBJ name to labels
    for e_ix = 1:numel(elec_sbj{sbj_ix}.label)
        elec_sbj{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix}.label{e_ix}];
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
    
    % Copy for other conditions
    for cond_ix = 2:numel(cond_lab)
        elec_sbj{sbj_ix,cond_ix} = elec_sbj{sbj_ix,1};
    end
    
    % Load Stats
    % Determine options: {'actv','CI','RT','CNI','pcon'}
    sig_ch = cell(size(cond_lab));
    if strcmp(stat_id,'actv')
        load([SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_actv_mn100.mat'],'actv_ch');
        sig_ch{1} = actv_ch;
        clear actv_ch
    elseif strcmp(stat_id,'CSE')
        load([SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'_',stat_id,'.mat'],'stat');
        for ch_ix = 1:numel(stat.label)
            if any(stat.mask(ch_ix,1,:))
                sig_ch{1} = [sig_ch{1} stat.label(ch_ix)];
            end
        end
        clear stat
    else    % ANOVA
        eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);        
        f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
        load(f_name,'stat','w2');
        
        % FDR correct pvalues for ANOVA
        for ch_ix = 1:numel(stat.label)
            pvals = squeeze(w2.pval(:,ch_ix,:));
            [~, ~, ~, qvals] = fdr_bh(pvals);%,0.05,'pdep','yes');
            
            % Consolidate to binary sig/non-sig
            for cond_ix = 1:numel(cond_lab)
                if strcmp(cond_lab{cond_ix},'RT') && any(stat.mask(ch_ix,1,:))
                    sig_ch{cond_ix} = [sig_ch{cond_ix} {[SBJs{sbj_ix} '_' stat.label{ch_ix}]}];
                elseif any(strcmp(cond_lab{cond_ix},{'CNI','pcon'})) && any(qvals(cond_ix,:)<0.05,2)
                    sig_ch{cond_ix} = [sig_ch{cond_ix} {[SBJs{sbj_ix} '_' w2.label{ch_ix}]}];
                end
            end
        end
        clear stat w2
    end
    
    % Remove hemi and/or atlas elecs
    if ~plot_out
        % Remove electrodes that aren't in atlas ROIs & hemisphere
        good_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, atlas_id, roi_id);
    else
        % Remove electrodes that aren't in hemisphere
        good_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, [], []);
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
            fprintf(sig_report,'%s - %s (%s)\n',sig_ch{cond_ix}{sig_ix},elec_sbj{sbj_ix,cond_ix}.atlas_lab{e_ix},...
                elec_sbj{sbj_ix,cond_ix}.hemi{e_ix});
        end
        
        % Select sig elecs && elecs matching atlas
        % fn_select_elec messes up if you try to toss all elecs
        good_elecs = intersect(good_elecs, sig_ch{cond_ix});
        if numel(intersect(elec_sbj{sbj_ix,cond_ix}.label,good_elecs))==0
            elec_sbj{sbj_ix,cond_ix} = {};
            good_sbj(sbj_ix,cond_ix) = false;
            warning('WARNING!!! All sig_ch are out of atlas and/or hemisphere!');
            fprintf(sig_report,'\t!!! 0 sig_ch remain\n');
        else
            cfgs = [];
            cfgs.channel = good_elecs;
            elec_sbj{sbj_ix,cond_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix,cond_ix});
            all_roi_labels{cond_ix} = [all_roi_labels{cond_ix}; elec_sbj{sbj_ix,cond_ix}.roi];
            all_roi_colors{cond_ix} = [all_roi_colors{cond_ix}; elec_sbj{sbj_ix,cond_ix}.roi_color];
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
    elec{cond_ix}.roi       = all_roi_labels{cond_ix};    % appendsens strips that field
    elec{cond_ix}.roi_color = all_roi_colors{cond_ix};    % appendsens strips that field
end

%% Load brain recon
mesh = fn_load_recon_mesh([],'mni',reg_type,hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
f = cell(size(cond_lab));
for cond_ix = 1%:numel(cond_lab)
    if strcmp(stat_id,'actv') || strcmp(stat_id,'CSE')
        plot_name = ['GRP_' stat_id '_' an_id];
    else
        plot_name = ['GRP_ANOVA_' cond_lab{cond_ix} '_' stat_id '_' an_id];
    end
    f{cond_ix} = figure('Name',plot_name);
    
    % Plot 3D mesh
    ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
    
    % Plot electrodes on top
    for e = 1:numel(elec{cond_ix}.label)
        cfgs = []; cfgs.channel = elec{cond_ix}.label(e);
        elec_tmp = fn_select_elec(cfgs,elec{cond_ix});
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.roi_color, 'label', lab_arg);
    end
    
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
    fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
        'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
        '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
    set(f{cond_ix}, 'windowkeypressfcn',   @cb_keyboard);
end


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
