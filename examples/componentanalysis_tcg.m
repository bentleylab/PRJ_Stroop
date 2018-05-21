function componentanalysis_tcg(subjectm, tcgver)

% Arjen's ICA computation script, 2018
% input should be 'subject_*' which is an mfile containing the info
% when using this script, you agree to supply him a lifetime's worth of
% red wine (Shiraz, Cab, or Pinot Noir)

% input check
if nargin == 0
  disp('Not enough input arguments');
  return;
end

% get subject info
eval(subjectm);

% preprocessing (continuous) data
data             = preproc_ica_tcg(subjectm, tcgver, true);

% downsample (still producing 60 Gb jobs)
if ~isfield(data, 'fsample') || ~isequal(data.fsample, 1000)
  cfg              = [];
  cfg.resamplefs   = 1000;
  cfg.detrend      = 'no';
  data             = ft_resampledata(cfg, data);
end

% component analysis
cfg              = [];
cfg.method       = 'runica';
comp             = ft_componentanalysis(cfg, data);
clear data

% store unmixing matrix and topolabel
comp             = keepfields(comp, {'unmixing', 'topolabel', 'topo', 'label'});
save([subj(tcgver).datadir subjectm(9:end) '_comp_' num2str(tcgver) '.mat'], 'comp');
clear comp
