function fn_view_recon_atlas_ROI(SBJ, pipeline_id, view_space, reg_type, show_labels, hemi, atlas_id, roi_id)
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

% View angle
if ~exist('view_angle','var')
    if strcmp(roi_id,'OFC')
        view_angle = [0 -90];   % from the bottom
    elseif strcmp(hemi,'l') && any([strcmp(roi_id,{'LPFC','INS','TMP','PAR','MTL'}) ~strcmp(roi_id,'MPFC')])
        view_angle = [-90 0];    % from the left
    elseif strcmp(hemi,'r') && any([strcmp(roi_id,{'LPFC','INS','TMP','PAR','MTL'}) ~strcmp(roi_id,'MPFC')])
        view_angle = [90 0];    % from the right
    else
        error(['Bad combo of hemi (' hemi ') and roi_id (' roi_id ')']);
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
    reg_suffix = ['_' reg_type];    % MNI space
else
    reg_suffix = '';                % Patient space
end

%% Load elec struct
[root_dir, ~] = fn_get_root_dir();
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

try
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_space,reg_suffix,'_',atlas_id,'_full.mat'];
    if exist([elec_fname(1:end-4) '_' roi_id '.mat'],'file')
        elec_fname = [elec_fname(1:end-4) '_' roi_id '.mat'];
    end
    load(elec_fname);
catch
    answer = input(['Could not load requested file: ' elec_fname ...
        '\nDo you want to run the atlas matching now? "y" or "n"\n'],'s');
    if strcmp(answer,'y')
        fn_save_elec_atlas(SBJ,pipeline_id,view_space,reg_type,atlas_id);
    else
        error('not running atlas assignment, exiting...');
    end
end

%% Match elecs to atlas ROIs
if any(strcmp(atlas_id,{'DK','Dx','Yeo7'}))
    if ~isfield(elec,'man_adj')
        elec.roi       = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
    end
    if strcmp(roi_id,'tissueC')
        elec.roi_color = fn_tissue2color(elec);
    elseif strcmp(atlas_id,'Yeo7')
        elec.roi_color = fn_atlas2color(atlas_id,elec.roi);
    else
        elec.roi_color = fn_roi2color(elec.roi);
    end
elseif any(strcmp(atlas_id,{'Yeo17'}))
    if ~isfield(elec,'man_adj')
        elec.roi       = elec.atlas_lab;
    end
    elec.roi_color = fn_atlas2color(atlas_id,elec.roi);
end

% Select relevant elecs (match ROI and hemi)
cfgs = [];
cfgs.channel = fn_select_elec_lab_match(elec, hemi, atlas_id, roi_id);
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
cfg.method      = 'iso2mesh';   % surface toolbox Arjen found
cfg.radbound    = 2;            % scalar indicating the radius of the target surface mesh element bounding sphere
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 100000;
cfg.smooth      = 3;
roi_mesh = ft_prepare_mesh(cfg, seg);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
h = figure;

% Plot 3D mesh
ft_plot_mesh(roi_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);

% Plot electrodes on top
cfgs = [];
for e = 1:numel(elec.label)
    cfgs.channel = elec.label(e);
    elec_tmp = fn_select_elec(cfgs, elec);
    ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.roi_color, 'label', lab_arg);
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(h, 'windowkeypressfcn',   @cb_keyboard);

%%
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
