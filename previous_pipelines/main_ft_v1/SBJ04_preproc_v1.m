function SBJ04_preproc(SBJ,line_filt)
%% Preprocess data using fieldtrip
% Inputs:
%   SBJ [str]- name of the subject
%   line_filt [str]- type of line noise filter
%       'bandstop'- traditional bandstop filter
%       'dft'- fieldtrip sinusoid regression
%       'cleanline'- Tim Mullen's cleanline function

% Parameters
psd_bp       = 'yes';          % plot psds after filtering?
psd_reref    = 'yes';          % plot psds after rereferencing?
psd_line     = 'yes';          % plot psds after filtering out line noise?
psd_fig_type = 'jpg';

demean_yn   = 'yes';
hp_yn       = 'yes';
lp_yn       = 'yes';
hp_freq     = 0.5;
hp_order    = 4;            % Leaving blank causes instability error, 1 or 2 works
lp_freq     = 300;

% line_filt   = 'bandstop';       % 'dft','bandstop', or 'cleanline'
% cleanline   = 'yes';                % Use Tim Mullen's cleanline function
% dft_yn      = 'no';
% bs_yn       = 'no';                % Parameters for this in preproc_vars

% Directories
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/'));
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

SBJ_dir     = ['/home/knight/hoycw/PRJ_Stroop/data/' SBJ '/'];
import_dir  = [SBJ_dir '01_import/'];
preproc_dir = [SBJ_dir '02_preproc/'];

%% Load data
fprintf('============== Loading Data %s ==============\n',SBJ);
import_filename = [import_dir SBJ '_1000hz.mat'];
load(import_filename);
data_orig = data;

preproc_vars_filename = [import_dir SBJ '_preproc_vars.mat'];
load(preproc_vars_filename);

%% High and low pass data
fprintf('============== High and Low Pass Filtering %s ==============\n',SBJ);
cfg           = [];
cfg.demean    = demean_yn;
cfg.lpfilter  = lp_yn;
cfg.lpfreq    = lp_freq;
cfg.hpfilter  = hp_yn;
cfg.hpfreq    = hp_freq;
cfg.hpfiltord = hp_order;
data = ft_preprocessing(cfg,data);
data_bp = data;

if strcmp(psd_bp,'yes')
    psd_dir = strcat(preproc_dir,'PSDs/bp/');
    if ~exist(psd_dir,'dir')
        mkdir(psd_dir);
    end
    fn_plot_PSD_1by1_compare_save(data_orig.trial{1},data.trial{1},data_orig.label,data.label,...
        data.fsample,strcat(psd_dir,SBJ,'_PSD_bp'),'orig','bp',psd_fig_type);
end
clear data_orig

%% Rerefrence
fprintf('============== Re-Referencing %s ==============\n',SBJ);
for d = 1:numel(preproc_vars.ch_lab.probes)
    cfg = [];
    cfg.channel = ft_channelselection(strcat(preproc_vars.ch_lab.probes{d},'*'), data.label);
    probe_data = ft_selectdata(cfg,data);   % Grab data from this probe to plot in PSD comparison
    
    cfg.montage.labelold = cfg.channel;
    [cfg.montage.labelnew, cfg.montage.tra, left_out_ch{d}] = fn_create_ref_scheme_bipolar(cfg.channel);
    data_reref{d} = ft_preprocessing(cfg, data);
    
    if strcmp(psd_reref,'yes')
        psd_dir = strcat(preproc_dir,'PSDs/bp.reref/');
        if ~exist(psd_dir,'dir')
            mkdir(psd_dir);
        end
        fn_plot_PSD_1by1_compare_save(probe_data.trial{1},data_reref{d}.trial{1},...
            probe_data.label,data_reref{d}.label,data_reref{d}.fsample,...
            strcat(psd_dir,SBJ,'_PSD_bp.reref'),'bp','bp.reref',psd_fig_type);
    end
end
clear data_bp probe_data

%Concatenate together again
cfg = [];
% cfg.appendsens = 'yes';
data = ft_appenddata(cfg,data_reref{:});
data_reref = data;

%% Filter out line noise
fprintf('============== Filtering Line Noise %s via %s ==============\n',SBJ,line_filt);
line_filt_id = strcat('bp.reref.',line_filt);

if strcmp(line_filt,'dft')
    cfg           = [];
    cfg.dftfilter = 'yes'; % line noise removal using discrete fourier transform
    cfg.dftfreq   = preproc_vars.notch_freqs;
    data = ft_preprocessing(cfg,data);
elseif strcmp(line_filt,'bandstop')
    % Calculate frequency limits
    bs_freq_lim = NaN([length(preproc_vars.notch_freqs) 2]);
    for f_ix = 1:length(preproc_vars.notch_freqs)
        bs_freq_lim(f_ix,:) = fn_freq_lim_from_CFBW(preproc_vars.notch_freqs(f_ix),preproc_vars.bs_width);
    end
    cfg          = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bs_freq_lim;
    data = ft_preprocessing(cfg,data);
elseif strcmp(line_filt,'cleanline')
    data = fn_cleanline(data,preproc_vars.notch_freqs);
else
    error('ERROR: line_filt type not in [dft, bandstop, cleanline]');
end

if strcmp(psd_line,'yes')
    psd_dir = strcat(preproc_dir,'PSDs/',line_filt_id,'/');
    if ~exist(psd_dir,'dir')
        mkdir(psd_dir);
    end
    fn_plot_PSD_1by1_compare_save(data_reref.trial{1},data.trial{1},data_reref.label,data.label,...
        data_reref.fsample,strcat(psd_dir,SBJ,'_PSD_',line_filt_id),'bp.reref',line_filt_id,psd_fig_type);
end
clear data_reref

%% Accumulate processing variables
proc_vars = preproc_vars;
proc_vars.ch_lab.left_out = [left_out_ch{:}];
proc_vars.line_filt_id    = line_filt_id;
proc_vars.line_type       = line_filt;
proc_vars.demean          = demean_yn;
proc_vars.hp              = hp_yn;
proc_vars.lp              = lp_yn;
proc_vars.hp_freq         = hp_freq;
proc_vars.hp_order        = hp_order;            % Leaving blank causes instability error, 1 or 2 works
proc_vars.lp_freq         = lp_freq;

%% Save data
output_filename = strcat(preproc_dir,SBJ,'_preproc_',line_filt_id,'.mat');
fprintf('============== Saving %s ==============\n',output_filename);
save(output_filename, '-v7.3', 'data');

proc_vars_filename = strcat(preproc_dir,SBJ,'_proc_vars.mat');
save(proc_vars_filename, '-v7.3', 'proc_vars');

end
