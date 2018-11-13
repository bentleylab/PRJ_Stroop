tick_step = 2;

cfgraw = [];
cfgraw.trials       = 'all';
%cfgraw.latency      = plt_lim;
cfgraw.channel      = 'all';
cfgraw.colorbar     = 'yes';
cfgraw.zlim         = 'maxabs'; %[-3 3];%
cfgraw.parameter    = 'powspctrm';
cfgraw.showlabels   = 'yes';
%cfgraw.baseline     = bsln_lim; %should be in sec
cfgraw.baselinetype = 'absolute';
%     cfg.layout = layout;
% cfg.maskparameter = 'trial';
% cfg.maskstyle = 'outline';

cfgdif = [];
cfgdif.trials       = 'all';
%cfgdif.channel      = 'all';
cfgdif.colorbar     = 'yes';
cfgdif.zlim         = 'maxabs'; %[-3 3];%
cfgdif.parameter    = 'powspctrm';
cfgdif.showlabels   = 'yes';
cfgdif.maskparameter = 'statmask';
cfgdif.maskstyle     = 'opacity';
cfgdif.maskalpha     = 0.5;
%cfgdif.latency       = stat_lim;
cfgdif.baseline      = [];
cfgdif.baselinetype  = [];
%cfgdif.channel       = tfr_diff.label{ch_ix};

