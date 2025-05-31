if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults


proc_id = 'main_ft';
view_space  = 'pat';
reg_type    = '';
atlas_id    = 'Dx';
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

%%
% Specify the filename
input_fname = '/Users/colinhoy/Code/PRJ_Stroop/Anas_Nicole Recon Changes.csv';

% Load the data from the CSV file into a table
opts = detectImportOptions(input_fname);      % Detect import options based on the file
opts = setvartype(opts, 'char');           % Set all variables to be of type 'char'
tbl = readtable(input_fname, opts);     % Read the file into a table

% Display the table
disp(tbl);
tbl.colin_roi = repmat({''},size(tbl.Label));
tbl.colin_groi = repmat({''},size(tbl.Label));
SBJs = unique(tbl.SubjectID);

%% Load elec and get ROI labels
for s = 1:length(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    elec_fname = [SBJ_vars.dirs.recon SBJs{s} '_elec_' proc_id '_' view_space reg_suffix '_orig_' atlas_id '_man.mat'];
    load(elec_fname);

    sbj_ix = find(strcmp(tbl.SubjectID,SBJs{s}));
    for ch_ix = 1:length(sbj_ix)
        tbl.colin_roi(sbj_ix(ch_ix)) = elec.ROI(strcmp(elec.label,tbl.Label(sbj_ix(ch_ix))));
        tbl.colin_groi(sbj_ix(ch_ix)) = elec.gROI(strcmp(elec.label,tbl.Label(sbj_ix(ch_ix))));
    end
end

tbl.groi_match = strcmp(tbl.NewGROI,tbl.colin_groi);

%% Write table with the manual labels added
output_fname = '/Users/colinhoy/Code/PRJ_Stroop/Anas_Nicole_Recon_Changes_CWH_edits.csv';
writetable(tbl, output_fname);

%% Compare Anas and my decisions
check_idx = ~strcmp(tbl.NewGROI,tbl.colin_groi) | ~strcmp(tbl.NewROI,tbl.colin_roi);
tbl_check = tbl(check_idx,:);

%% Plot electrodes to check
SBJ = 'IR413';
probe_prefix = 'LAC';

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix '_orig_' atlas_id '_man.mat'];
load(elec_fname);
elec.color = fn_roi2color(elec.ROI);
pial_l = fn_load_recon_mesh(SBJ,view_space,reg_type,'pial','l');
pial_r = fn_load_recon_mesh(SBJ,view_space,reg_type,'pial','r');
wm_l   = fn_load_recon_mesh(SBJ,view_space,reg_type,'wm','l');
wm_r   = fn_load_recon_mesh(SBJ,view_space,reg_type,'wm','r');

cfgs.channel = ft_channelselection([probe_prefix '*'], elec.label);
probe = fn_select_elec(cfgs, elec);

% Plot pial
mesh_alpha = 0.7;
if numel(probe.label)>1; hemi = probe.hemi{2}; else hemi = probe.hemi{1};end
pial = figure('Name',[SBJ ' pial ' probe.label{1} ' : ' probe.label{end}],...
    'Units','normalized');
set(pial,'OuterPosition',[0 0 0.5 1]);
ft_plot_mesh(eval(['pial_' hemi]), 'vertexcolor', 'curv', 'facealpha', mesh_alpha);
for e = 1:numel(probe.label)
    cfgs.channel = probe.label(e);
    tmp = fn_select_elec(cfgs, probe);
    ft_plot_sens(tmp, 'elecshape', 'sphere', 'facecolor', tmp.color, 'label', 'label');
end
material dull; lighting gouraud;
l = camlight;
set(pial, 'windowkeypressfcn',   @cb_keyboard);

% Plot WM
if strcmp(SBJ_vars.ch_lab.probe_type,'seeg')
    wm = figure('Name',[SBJ ' white ' probe.label{1} ' : ' probe.label{end}],...
        'Units','normalized');
    set(wm, 'OuterPosition', [0.5 0 0.5 1]);
    ft_plot_mesh(eval(['wm_' hemi]), 'vertexcolor', 'curv', 'facealpha', 0.5);
    for e = 1:numel(probe.label)
        cfgs.channel = probe.label(e);
        tmp = fn_select_elec(cfgs, probe);
        ft_plot_sens(tmp, 'elecshape', 'sphere', 'facecolor', tmp.color, 'label', 'label');
    end
    material dull; lighting gouraud;
    l = camlight;
    set(wm, 'windowkeypressfcn',   @cb_keyboard);
end

% fn_view_recon(SBJ, '', 'ortho', view_space, reg_type, 1, 'b', 1);
