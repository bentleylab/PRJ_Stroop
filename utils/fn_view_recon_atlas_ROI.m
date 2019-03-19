function fn_view_recon_atlas_ROI(SBJ, pipeline_id, view_space, reg_type, show_labels, hemi, atlas_id, roi_id)%, view_angle)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJ [str] - subject ID to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   plot_type [str] - {'ortho', '3d'} choose 3 slice orthogonal plot or 3D surface rendering
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r'} hemisphere to plot (can't be both, that's a shitty plot)
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - gROI grouping to pick mesh and color specific ROIs
%       'LPFC','MPFC','OFC','INS','TMP','PAR'

%% Process Inputs
% Error cases
if strcmp(hemi,'b') && ~strcmp(roi_id,'OFC')
    error('hemi must be l or r for all non-OFC plots');
end
if ~any(strcmp(roi_id,{'LPFC','MPFC','INS','OFC','TMP','PAR'}))
    error('roi_id needs to be a lobe (not OCC either)');
end

% View angle
if strcmp(roi_id,'OFC')
    view_angle = [0 -90];   % from the bottom
elseif strcmp(hemi,'l') && any([strcmp(roi_id,{'LPFC','INS','TMP','PAR','MTL'}) ~strcmp(roi_id,'MPFC')])
    view_angle = [-90 0];    % from the left
elseif strcmp(hemi,'r') && any([strcmp(roi_id,{'LPFC','INS','TMP','PAR','MTL'}) ~strcmp(roi_id,'MPFC')])
    view_angle = [90 0];    % from the right
else
    error(['Bad combo of hemi (' hemi ') and roi_id (' roi_id ')']);
%     view_angle = [0 0]; %straight on
end

% Suffixes
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end
if strcmp(roi_id,'tissue') || strcmp(roi_id,'tissueC')
    tis_suffix = '_tis';
else
    tis_suffix = '';
end

%% Load elec struct
[root_dir, ~] = fn_get_root_dir();
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

try
    elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_space,reg_suffix,'_',atlas_id,tis_suffix,'.mat'];
    load(elec_atlas_fname);
catch
    answer = input(['Could not load requested file: ' elec_atlas_fname ...
        '\nDo you want to run the atlas matching now? "y" or "n"\n'],'s');
    if strcmp(answer,'y')
        fn_save_elec_atlas(SBJ,pipeline_id,view_space,reg_type,atlas_id);
    else
        error('not running atlas assignment, exiting...');
    end
end

%% Match elecs to atlas ROIs
[roi_list, ~] = fn_roi_label_styles(roi_id);

if any(strcmp(atlas_id,{'DK','Dx','Yeo7'}))
    elec.roi       = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
    if strcmp(roi_id,'tissueC')
        elec.roi_color = fn_tissue2color(elec);
    elseif strcmp(atlas_id,'Yeo7')
        elec.roi_color = fn_atlas2color(atlas_id,elec.roi);
    else
        elec.roi_color = fn_roi2color(elec.roi);
    end
elseif any(strcmp(atlas_id,{'Yeo17'}))
    elec.roi       = elec.atlas_lab;
    elec.roi_color = fn_atlas2color(atlas_id,elec.roi);
end

% Find elecs matching ROI
roi_match = false([numel(elec.label) numel(roi_list)]);
for roi_ix = 1:numel(roi_list)
    roi_match(:,roi_ix) = strcmp(elec.roi,roi_list{roi_ix});
end
% Exclude other hemisphere
if ~strcmp(hemi,'b')
    hemi_out_elecs = elec.label(~strcmp(elec.hemi,hemi));
else
    hemi_out_elecs = {};
end

% Select relevant elecs
cfgs = [];
cfgs.channel = [elec.label(any(roi_match,2)); fn_ch_lab_negate(hemi_out_elecs)'];
elec = fn_select_elec(cfgs,elec);

%% Load Atlas
atlas = fn_load_recon_atlas(SBJ,atlas_id);

% Get Atlas-ROI mapping
atlas_labels = fn_atlas_roi_select_mesh(atlas_id, roi_id, hemi);

%% Select ROI mesh
cfg = [];
cfg.inputcoord = atlas.coordsys;
cfg.atlas = atlas;
cfg.roi = atlas_labels;
roi_mask = ft_volumelookup(cfg,atlas);

seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = roi_mask;

cfg = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;            % scalar indicating the radius of the target surface mesh element bounding sphere
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 100000;
cfg.smooth      = 3;
roi_mesh = ft_prepare_mesh(cfg, seg);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
h = figure;

% Plot 3D mesh
mesh_alpha = 0.8;
if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
    mesh_alpha = 0.3;
end
ft_plot_mesh(roi_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);

% Plot electrodes on top
cfgs = [];
for e = 1:numel(elec.label)
    cfgs.channel = elec.label(e);
    elec_tmp = fn_select_elec(cfgs, elec);
    if show_labels
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.roi_color, 'label', 'label');
    else
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.roi_color);
    end
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(h, 'windowkeypressfcn',   @cb_keyboard);

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
