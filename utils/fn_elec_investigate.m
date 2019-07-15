function fn_elec_investigate()

% Load elec

% Load mri

% Load meshes
% pial
% smoothwm
% inflated + coloring

% Plot meshes
i = figure('Name','inflated');
ft_plot_mesh(inflated, 'vertexcolor', 'curv');
material dull; lighting gouraud;
l = camlight;
set(i, 'windowkeypressfcn',   @cb_keyboard);

end