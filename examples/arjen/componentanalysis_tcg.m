function componentanalysis_tcg(subjectm, role)
% 'memreq', 15*1024^3, 'timreq', 8*3600

% this is an ICA script
% input should be 'subject_pair*' which is an mfile containing the info
% role should be 'sender' or 'receiver'
% A. Stolk

% input check
if nargin == 0
  disp('Not enough input arguments');
  return;
end

% get subject info
eval(subjectm);

% preprocessing
data = preproc_tcg(subjectm, true); % use all/continuous data for ICA

% select sender/receiver EEG channels
if strcmp(role, 'sender')
  cfg              = [];
  cfg.channel      = ft_channelselection({'1-F*','1-A*','1-C*','1-T*','1-P*','1-O*','1-I*'}, data.label);
  eeg              = ft_selectdata(cfg, data);
  rank             = numel(eeg.label)-1; % minus common average
elseif strcmp(role, 'receiver')
  cfg              = [];
  cfg.channel      = ft_channelselection({'2-F*','2-A*','2-C*','2-T*','2-P*','2-O*','2-I*'}, data.label);
  eeg              = ft_selectdata(cfg, data);
  rank             = numel(eeg.label)-1; % minus common average
end
clear data

% component analysis
cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = rank;
ica              = ft_componentanalysis(cfg, eeg);
clear eeg

% keep unmixing matrix and topolabel
comp.unmixing    = ica.unmixing;
comp.topolabel   = ica.topolabel;
comp.topo        = ica.topo;
comp.label       = ica.label;
clear ica

% store output
save([subj(1).datadir subjectm(9:end) '_' role(1) 'comp.mat'], 'comp'); % scomp or rcomp
clear comp