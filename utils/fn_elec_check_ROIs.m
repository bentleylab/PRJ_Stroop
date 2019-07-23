function fn_elec_check_ROIs(SBJ)%, proc_id, view_space, reg_type, atlas_id)
%% Plot recons for each probe with pial, white matter, and inflated to check ROI assignments
proc_id = 'main_ft';
view_space  = 'pat';
reg_type    = '';
atlas_id    = 'Dx';

%% Load elec and get ROI labels
% Check which root directory
[root_dir, ~] = fn_get_root_dir();
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space '_' reg_suffix atlas_id '_compiled.mat'];
load(elec_fname);

% Assign ROIs
elec.color = fn_roi2color(elec.ROI);

%% Load meshes
pial_l = fn_load_recon_mesh(SBJ,view_space,reg_type,'pial','l');
pial_r = fn_load_recon_mesh(SBJ,view_space,reg_type,'pial','r');
wm_l   = fn_load_recon_mesh(SBJ,view_space,reg_type,'wm','l');
wm_r   = fn_load_recon_mesh(SBJ,view_space,reg_type,'wm','r');
infl_l = fn_load_recon_mesh(SBJ,view_space,reg_type,'inflated','l');
infl_r = fn_load_recon_mesh(SBJ,view_space,reg_type,'inflated','r');

% Plot inflated (no elecs)
i_l = figure('Name',[SBJ ' inflated L']);
ft_plot_mesh(infl_l, 'vertexcolor', 'curv');
material dull; lighting gouraud;
l = camlight;
set(i_l, 'windowkeypressfcn',   @cb_keyboard);

i_r = figure('Name',[SBJ ' inflated L']);
ft_plot_mesh(infl_r, 'vertexcolor', 'curv');
material dull; lighting gouraud;
l = camlight;
set(i_r, 'windowkeypressfcn',   @cb_keyboard);

%% Plot elecs + meshes
for p = 1:numel(SBJ_vars.ch_lab.probes)
    cfgs.channel = ft_channelselection([SBJ_vars.ch_lab.probes{p} '*'], elec.label);
    probe = fn_select_elec(cfgs, elec);
    
    if strcmp(SBJ_vars.ch_lab.probe_type{p},'ecog')
        mesh_alpha = 0.8;
    else
        mesh_alpha = 0.3;
    end
    
    % Plot pial
    pial = figure('Name',[SBJ ' pial ' probe.label{1} ' : ' probe.label{end}]);
    ft_plot_mesh(eval(['pial_' probe.hemi{2}]), 'vertexcolor', 'curv', 'facealpha', mesh_alpha);
    for e = 1:numel(probe.label)
        cfgs.channel = probe.label(e);
        tmp = fn_select_elec(cfgs, probe);
        ft_plot_sens(tmp, 'elecshape', 'sphere', 'facecolor', tmp.color);
    end
    material dull; lighting gouraud;
    l = camlight;
    set(pial, 'windowkeypressfcn',   @cb_keyboard);
    
    % Plot WM
    if strcmp(SBJ_vars.ch_lab.probe_type,'seeg')
        wm = figure('Name',[SBJ ' white ' probe.hemi{2}]);
        ft_plot_mesh(eval(['wm_' probe.hemi{2}]), 'vertexcolor', 'curv', 'facealpha', 0.5);
        for e = 1:numel(probe.label)
            cfgs.channel = probe.label(e);
            tmp = fn_select_elec(cfgs, probe);
            ft_plot_sens(tmp, 'elecshape', 'sphere', 'facecolor', tmp.color);
        end
        material dull; lighting gouraud;
        l = camlight;
        set(wm, 'windowkeypressfcn',   @cb_keyboard);
    end
end

%% Plot Elecs
if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
    fn_view_recon(SBJ, proc_id, 'ortho', view_space, reg_type, 1, 'b', 1);
end

end