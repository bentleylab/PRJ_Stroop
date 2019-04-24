if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Run fn_view_recon_stat_movie
% General parameters
SBJs = {'CP24','CP26','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
        'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};
pipeline_id = 'main_ft';
plt_id      = 'recon_movie';
reg_type    = '';
view_space  = 'pat';
plot_out = 0;
hemi = 'r';

save_fig    = 1;
fig_vis     = 'on';
fig_type    = 'svg';

% Analysis parameters
actv_win    = '100';



stat_id = 'corrRT_CNI_pcon_WL200_WS50';

%% R-locked
an_id = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';

for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fn_view_recon_stat_movie(SBJ, pipeline_id, stat_id, an_id, view_space, reg_type, 'r', 0, plt_id);
%     fn_view_recon_stat_movie(SBJ, pipeline_id, stat_id, an_id, view_space, reg_type, 'l', 0, plt_id, varargin);
end
