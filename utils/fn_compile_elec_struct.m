function fn_compile_elec_struct(SBJ,pipeline_id,view_space)
%% Compile ROI info from single electrodes into bipolar pairs
%   Goes from smallest to largest (ELEC1-ELEC2, ELEC2-ELEC3, etc.)
%   Pairs are drawn from preprocessed data labels
%
% new_labels = cell array of strings with combined labels
%   can be used with ft_preprocessing as cfg_rereference.montage.labelnew
% weights = [length(new_labels),length(labels)] matrix of weights to combine elecs
%   can be used with ft_preprocessing as cfg_rereference.montage.tra
% left_out_ch = labels of channels that aren't included (don't have contiguous pair)

% Set up paths
[root_dir, ft_dir] = fn_get_root_dir();
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load data
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{1});
else
    block_suffix = SBJ_vars.block_name{1};   % should just be ''
end
import_filename = [SBJ_vars.dirs.import SBJ '_',num2str(proc_vars.resample_freq),'hz',block_suffix,'.mat'];
load(import_filename);

% % Original (single electrode) labels
% import  = load([SBJ_vars.dirs.import SBJ '_1000hz.mat']);
% raw_lab = import.data.label;

%% Load Elec struct
if strcmp(view_space,'pat')
    elec_name = SBJ_vars.recon.elec_pat;
elseif strcmp(view_space,'mni')
    elec_name = SBJ_vars.recon.elec_mni;
else
    error(['Unknown view_space: ' view_space]);
end
slash = strfind(elec_name,'/'); elec_suffix = elec_name(slash(end)+numel(SBJ)+2:end-4);
tmp = load(elec_name); eval(['elec = tmp.' elec_suffix ';']); clear tmp;

%% Select imported channels
cfg = []; cfg.channel = data.label;
elec = fn_select_elec(cfg,elec);

%% Apply montage per probe
left_out_ch = {};
elec_reref  = cell([1 numel(SBJ_vars.ch_lab.probes)]);
for d = 1:numel(SBJ_vars.ch_lab.probes)
    cfg = [];
    cfg.channel = ft_channelselection(strcat(SBJ_vars.ch_lab.probes{d},'*'), data.label);
    probe_elec  = fn_select_elec(cfg,elec);
    %     probe_data = ft_selectdata(cfg,data);   % Grab data from this probe to plot in PSD comparison
    %     probe_data.elec = fn_elec_ch_select(elec,cfg.channel);
    
    % Create referencing scheme
    if strcmp(SBJ_vars.ref_types{d},'BP')
        cfg.montage.labelold = cfg.channel;
        [cfg.montage.labelnew, cfg.montage.tra, left_out_ch{d}] = fn_create_ref_scheme_bipolar(cfg.channel);
        cfg.updatesens = 'yes';
        elec_reref{d} = ft_apply_montage(probe_elec, cfg.montage);%, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
        %     data_reref{d} = ft_preprocessing(cfg, probe_data);
    else
        elec_reref{d} = probe_elec;
    end
end

%% Recombine
cfg = [];
elec = ft_appendsens(cfg,elec_reref{:});
elec.type = 'ieeg';
for e = 1:numel(elec.chantype)
    elec.chantype{e} = lower(SBJ_vars.ch_lab.probe_type{d});
end

%% Save data
output_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_elec_',pipeline_id,'_',view_space,'.mat');
fprintf('============== Saving %s ==============\n',output_filename);
save(output_filename, '-v7.3', 'elec');

end
