function componentdetection_tcg(subjectm, role)

% this is an interactive preprocessing+cleaning script
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
inputfile = [subj(1).datadir subjectm(9:end) '_' role(1) 'comp.mat'];

% check for prior detection and whether the unmixing matrix is available
if ~isfield(subj(1), ['bad' role(1) 'comp']) && ~isempty(dir(inputfile))
  
  % preprocessing
  [data, ext] = preproc_tcg(subjectm); % trial-segmented data
  
  % select sender/receiver EEG and EOG channels
  if strcmp(role, 'sender')
    cfg              = [];
    cfg.channel      = ft_channelselection({'1-F*','1-A*','1-C*','1-T*','1-P*','1-O*','1-I*'}, data.label);
    eeg              = ft_selectdata(cfg, data);
    cfg.channel      = {'1-EXG1','1-EXG2','1-EXG3','1-EXG4'}; % above, below, aside right eye, and aside left eye
    ext              = ft_selectdata(cfg, ext);
  elseif strcmp(role, 'receiver')
    cfg              = [];
    cfg.channel      = ft_channelselection({'2-F*','2-A*','2-C*','2-T*','2-P*','2-O*','2-I*'}, data.label);
    eeg              = ft_selectdata(cfg, data);
    cfg.channel      = {'2-EXG1','2-EXG2','2-EXG3','2-EXG4'}; % above, below, aside right eye, and aside left eye
    ext              = ft_selectdata(cfg, ext);
  end
  for t = 1:numel(ext.trial)
    eog.trial{1,t}(1,:) = ext.trial{1,t}(1,:) - ext.trial{1,t}(2,:); % vertical EOG trace
    eog.trial{1,t}(2,:) = ext.trial{1,t}(3,:) - ext.trial{1,t}(4,:); % horizontal EOG trace
  end
  clear data ext

  % decompose the datasegments into components, using the previously found (un)mixing matrix
  load(inputfile, 'comp');
  cfg           = [];
  cfg.unmixing  = comp.unmixing;
  cfg.topolabel = comp.topolabel;
  ica           = ft_componentanalysis(cfg, eeg);
  clear comp eeg
  
  % calculate correlation with EOG traces
  eogcorr = eogcorrelation(ica, eog);
  display(eogcorr);
  clear eog eogcorr
  
  % calculate powerspectra
  cfg           = [];
  cfg.method    = 'mtmfft';
  cfg.taper     = 'hanning';
  cfg.foilim    = [1 200];
  freq          = ft_freqanalysis(cfg, ica);
  for c = 1:ceil(numel(freq.label)/8) % per 8 units
    cmps = c*8-7:c*8;
    while cmps(end) > numel(freq.label) % in case fewer channels than 64
      cmps(end) = [];
    end
    figure; semilogy(freq.freq, mean(freq.powspctrm,1), 'k', 'linewidth', 2)
    hold on; semilogy(freq.freq, freq.powspctrm(cmps,:), 'linewidth', 1)
    legend('avg', freq.label{cmps})
    xlim([1 200])
  end
  
  % visualize topographies and signals
  if strcmp(role, 'sender')
    ica.topolabel   = regexprep(ica.topolabel, '1-', ''); % remove '1-' prefix
  elseif strcmp(role, 'receiver')
    ica.topolabel   = regexprep(ica.topolabel, '2-', ''); % remove '2-' prefix
  end
  cfg           = [];
  cfg.viewmode  = 'component';
  cfg.layout    = 'biosemi64.lay'; % specify the layout file that should be used for plotting
  cfg.fontsize  = 1/numel(ica.topolabel);
  cfg.linewidth = 1;
  ft_databrowser(cfg, ica);
  clear ica
  
  % note bad components
  rejcomp = input('specify components for rejection (e.g. [1 2 6]): ');
  close all
  
  % enter visually detected artifacts in subject m-file
  fid = fopen([subjectm '.m'],'At');
  fprintf(fid,'\n%s\n',['%% Entered @ ' datestr(now)]);
  fprintf(fid,'%s\n', ['subj(1).bad' role(1) 'comp = [' num2str(rejcomp) '];']);
  fclose all;
  
end


%%% SUBFUNCTION %%%
function [eogcorr] = eogcorrelation(comp, eog)

% comp               - componentanalyzed EEG/MEG chans
% eog                - preprocessed EOG chans
% Arjen Stolk

Ncomp = numel(comp.label);
eogcorr = cat(2,[1:1:Ncomp]',zeros(Ncomp,1));
Ntrls = length(comp.trial);

for i=1:Ntrls
    
    % compute orthonormal basis for subspace spanned by artefact signals
    X = eog.trial{i}'; % note the transpose
    % remove the mean
    X = X - ones(size(X,1),1)*mean(X,1);
    % compute Singular value decomposition (svd) of X
    [U,S,V] = svd(X,0);
    % the left singular vectors ('U') contain the orthonormal basis
       
    % compute subspace correlation of each component waveform with U
    for j=1:Ncomp
        % get component wavefrom and standardize
        x1 = comp.trial{i}(j,:);
        x1 = x1(:);
        x1 = x1 - mean(x1);
        x1 = x1./norm(x1,2);
        % correlate
        eogcorr(j,2) = eogcorr(j,2) + sqrt(sum(abs(U'*x1).^2));
    end
end

% divide by number of trials, and sort on basis of correlation
eogcorr(:,2) = eogcorr(:,2)./Ntrls;
eogcorr = sortrows(eogcorr,-2);