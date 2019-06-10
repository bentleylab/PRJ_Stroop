plt_vars.fig_pos  = {[0 0 0.5 0.7], [0 0 0.5 0.8]};
plt_vars.axis_vis = 'off';

% Venn features
plt_vars.cmap       = parula(100);
plt_vars.ang        = linspace(0, 2*pi);
plt_vars.ax_fudge   = 3;
plt_vars.alpha      = 0.6;
plt_vars.line_color = [0 0 0];

% Text features
plt_vars.title_sz     = 20;
plt_vars.legend_sz    = 16;
plt_vars.text_sz      = 20;
plt_vars.text_align   = 'center';
plt_vars.text_spacers = {...
        [-1 0; 1 0; ...%Main Effects
         0 -1],  % single pair
        [-1 0; 1 0; 0 1; ...%Main Effects
         0 -1; -1 1; 1 1; 0 0]};  % pairs + interaction
