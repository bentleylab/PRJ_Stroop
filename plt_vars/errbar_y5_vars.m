plt.fig_pos    = [0 0 0.5 0.6];
plt.legend_loc = 'northeast';

plt.errbar_color = 'k';
plt.errbar_width = 1.5;
plt.scat_style   = 'k*';
plt.scat_sz      = 50;

plt.scat_offsets = linspace(-0.015,0.015,numel(SBJs));
plt.bar_offsets  = linspace(-0.25,0.25,numel(roi_list));

plt.ax_color    = 'k';
plt.ylim        = [0 0.5];
plt.y_tick_step = 0.1;
plt.ylab_sz     = 16;
plt.title_sz    = 20;
