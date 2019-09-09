plt_vars.movie_lim  = [-1 -1];  % should take -50ms to end
plt_vars.play_speed = 0.5;

% Electrode properties
plt_vars.cmap_name = 'parula';      % colormap for electrodes
plt_vars.ns_color  = [0 0 0];       % color of non-significant spheres (black)
plt_vars.ns_sz     = 0.5;           % radius of non-significant spheres
plt_vars.sz_lim    = [1 5];         % radius (in mm) of spheres

% Event + Time lines
plt_vars.evnt_type          = {'S'};
plt_vars.evnt_len           = 0.1;  % time to reamin on screen
plt_vars.evnt_str_pos       = [0 0 -110]; % X, Y, Z
plt_vars.evnt_str_color     = 'k';
plt_vars.evnt_str_size      = 12;
plt_vars.evnt_str_horzalign = 'center';
plt_vars.evnt_str_weight    = 'bold';
%plt_vars.time_ln_color = 'k';
%plt_vars.time_ln_width = 1;
%plt_vars.time_ln_len   = 150;
%plt_vars.time_ln_z     = -110;

% Time String
plt_vars.time_str_pos       = [0 0 -100]; % X, Y, Z
plt_vars.time_str_color     = 'k';
plt_vars.time_str_size      = 12;
plt_vars.time_str_horzalign = 'center';

% Movie properties
plt_vars.vid_ext      = '.mp4';
plt_vars.vid_encoding = 'MPEG-4';
plt_vars.frame_delay  = 0.02;       % delay between saving frames (allows for monitor refresh)
plt_vars.frame_skip   = 1;          % 1 = plot every time point, 2 = every other, etc.
