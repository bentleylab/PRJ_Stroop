%% Regina 3D mesh and electrode plotting example code

%% Load elec
elec_fname = 'your_elec_file.mat';
load(elec_fname);
% maybe also add an elec.color, .roi, or some such variable with the automatic
% (maximum probability) assignments, which can be used to color the elecs
% (see elec_rgb below) and then double check the accuracy visually

%% Load atlas/ROIs
roi_mr = ft_read_mri('file.nii');
% This should have an [x,y,x] sized field of values that are 0/1
%   Not sure if you need to, but I would convert it to logicals
%   note the name of this field for .tissue below, likely roi_mr.anatomy
% make sure it has a .coordsys field
roi_mr.coordsys = 'acpc';   % this should work, and doesn't matter too much for this case

%% Create mesh of ROI
cfg = [];
cfg.method      = 'iso2mesh';   % surface toolbox Arjen found
cfg.radbound    = 2;            % scalar indicating the radius of the target surface mesh element bounding sphere
cfg.maxsurf     = 0;
cfg.tissue      = 'anatomy';    % this should match the fieldname of your [x,y,z] logical field
cfg.numvertices = 100000;       % increase for nicer (slower) meshes; 100k is probably overkill for just amygdala
cfg.smooth      = 3;
roi_mesh = ft_prepare_mesh(cfg, roi_mr);

%% Plot mesh and electrodes
h = figure;

% Plot 3D mesh
mesh_alpha = 0.3;   % decent for SEEG
ft_plot_mesh(roi_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);

% Plot all electrodes with single color
elec_rgb = [1 0 0];  %make this whatever you want
lab_arg = 'label';      % or 'off' for no labels
ft_plot_sens(elec, 'elecshape', 'sphere', 'facecolor', elec_rgb, 'label', lab_arg);

% % Alternative way to plot elecs:
% % Get X,Y,Z of sphere points
% [xsp, ysp, zsp] = sphere(100);
% % Plot as a surface (and scale to size)
% size = 1; % in mm (I think) if you are plotting it in the same coordsys as an MR surface mesh
% for e_ix = 1:numel(elec.label)
%     e_sphr{e_ix} = surf(size*xsp+elec.chanpos(e_ix,1), size*ysp+elec.chanpos(e,2), size*zsp+elec.chanpos(e,3));
%     set(e_sphr{e_ix}, 'EdgeColor', elec_rgb);
% end

% Make the lighting nice
%   (you might need to paste the 2 functions below into their own .m files)
view(view_angle); material dull; lighting gouraud;
% enable option to press "l" to make the camera light shine
%   from the current view angle
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(h, 'windowkeypressfcn',   @cb_keyboard);


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

