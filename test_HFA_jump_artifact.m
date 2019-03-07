%% test filtering interactions for the jump in HFA
SBJ = 'IR31';
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
import = load([SBJ_vars.dirs.import SBJ '_1000hz.mat']);
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

cfgs = []; cfgs.channel = {'ROF1','ROF2'};
data = ft_selectdata(cfgs,import.data);

cfg = [];
% cfg.lpfilter = 'yes';
% cfg.lpfreq   = 300;
% cfg.hpfilter = 'yes';
% cfg.hpfreq   = 0.5;
% cfg.hpfiltord = 6;
cfg.bpfilter = 'yes';
cfg.bpfreq = [0.5 300];
cfg.demean   = 'yes';
data = ft_preprocessing(cfg,data);

data_both = data;
cfgs.channel = {'ROF1'};
data = ft_selectdata(cfgs,data);
data.trial{1} = data.trial{1}-data_both.trial{1}(2,:);
data.label{1} = 'ROF1-2';

bs_freq_lim = NaN([length(SBJ_vars.notch_freqs) 2]);
for f_ix = 1:length(SBJ_vars.notch_freqs)
    bs_freq_lim(f_ix,:) = fn_freq_lim_from_CFBW(SBJ_vars.notch_freqs(f_ix),SBJ_vars.bs_width);
end
bs_freq_lim(bs_freq_lim(:,2) > data.fsample/2, :) = [];
cfg          = [];
cfg.bsfilter = 'yes';
cfg.bsfreq   = bs_freq_lim;
data = ft_preprocessing(cfg,data);

ft_databrowser(cfgp,data);

%%
roi = data;
trial_lim_s_pad = [-1 3];
bsln_events = trial_info.word_onset;
max_RT  = max(trial_info.response_time);
roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,trial_info.condition_n',...
    round(trial_lim_s_pad*data.fsample));

%%
HFA_type   = 'hilbert';
foi_lim = [70 150]; % min and max of desired frequencies
n_foi   = 8;
min_exp = log(foi_lim(1))/log(2); % that returns the exponents
max_exp = log(foi_lim(2))/log(2);
fois    = 2.^[linspace(min_exp,max_exp,n_foi)];
foi_bws = fn_semilog_bws(fois);     % semilog bandwidth spacing to match Erik Edwards & Chang lab
bp_lim  = zeros([numel(fois) 2]);
for f = 1:numel(fois)
    bp_lim(f,:) = fn_freq_lim_from_CFBW(fois(f), foi_bws(f));
end

cfg_hfa = [];
cfg_hfa.hilbert  = 'abs';
cfg_hfa.bpfilter = 'yes';
cfg_hfa.bpfreq   = [];      % to be filled by looping through foi_center
cfg_hfa.channel  = 'all';

hfas = cell(size(fois));
orig_lab = roi_trl.label;
for f_ix = 1:numel(fois)
    cfg_hfa.bpfreq = bp_lim(f_ix,:);
    cfg_hfa.demean = 'yes';
    fprintf('\n------> %s filtering: %.03f - %.03f\n', HFA_type, bp_lim(f_ix,1), bp_lim(f_ix,2));
    hfas{f_ix} = ft_preprocessing(cfg_hfa,roi_trl);
    hfas{f_ix}.label = strcat(hfas{f_ix}.label,[':' num2str(fois(f_ix),4)]);
end
% Treat different freqs as channels
hfa_tmp = hfas{1};      % Save to plug in averaged data
hfa = ft_appenddata([], hfas{:}); clear hfas;

%% 
cfg_class = [];
cfg_class.bpfilter = 'yes';
cfg_class.bpfreq = [70 150];
cfg_class.hilbert = 'abs';
class = ft_preprocessing(cfg_class,roi_trl);

%%
cfg_plot = [];
cfg_plot.viewmode = 'vertical';
ft_databrowser(cfg_plot,hfa);

%%
cfgd = [];
cfgd.derivative = 'yes';
hfa_d = ft_preprocessing(cfgd,hfa);
ft_databrowser(cfg_plot,hfa_d);

 %%
figure;
for ch = 1:numel(hfa_d.label)
    plot(hfa.trial{1}(ch,:));
    pause;
end