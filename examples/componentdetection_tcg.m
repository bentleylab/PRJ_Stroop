function componentdetection_tcg(subjectm, tcgver)

% Arjen's ICA detection script, 2018
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
inputfile = [subj(tcgver).datadir subjectm(9:end) '_comp_' num2str(tcgver) '.mat'];

% check for prior detection and whether the unmixing matrix is available
if ~isfield(subj(tcgver), 'badcomp') && ~isempty(dir(inputfile))
  
  % create output folder
  if ~exist([subj(tcgver).datadir 'ica'], 'dir')
    mkdir([subj(tcgver).datadir 'ica']);
  end
  
  % use the previously found (un)mixing matrix
  load(inputfile, 'comp');
  
  % component uniformity
  for c = 1:size(comp.topo,2)
    % all elec
    observed         = abs(comp.topo(:,c)); % channel weights are the counts
    expected         = mean(observed); % uniform weight / count
    comp.stat(c,1)   = sum((observed - expected).^2 / expected); % pearson chi-squared value
    comp.prob(c,1)   = chi2cdf(comp.stat(c,1), size(comp.topo,1)-1, 'upper'); % p-value (larger ~ spatially more broad)
    comp.critval     = chi2inv(1-0.05, size(comp.topo,1)-1); % critval for the dofs 
  end
  
  % get electrode information
  if isequal(exist([subj(tcgver).datadir 'recon/' subjectm(9:end) '_elec_acpc_fr.mat'], 'file'),2)
    load([subj(tcgver).datadir 'recon/' subjectm(9:end) '_elec_acpc_fr.mat']); % fused & realigned
    elec = elec_acpc_fr; clear elec_acpc_fr
  else
    load([subj(tcgver).datadir 'recon/' subjectm(9:end) '_elec_acpc_f.mat']); % fused
    elec = elec_acpc_f; clear elec_acpc_f
  end
  
  % create topography
  mesh           = ft_read_headshape([subj(tcgver).datadir '/recon/freesurfer/surf/lh.pial']);
  mesh.coordsys  = 'acpc';
  comp.label     = comp.topolabel; % needed for layout creation
  cfg            = [];
  cfg.elec       = elec;
  cfg.headshape  = mesh;
  cfg.projection = 'orthographic';
  cfg.channel    = 'all';
  cfg.viewpoint  = 'left';
  cfg.mask       = 'convex';
  lay            = ft_prepare_layout(cfg, comp); % ft_plot_lay(lay)
  
  % preprocessing
  data           = preproc_ica_tcg(subjectm, tcgver, false); % trial-based data
  
  % decompose the datasegments into components, using the previously found unmixing matrix
  cfg           = [];
  cfg.unmixing  = comp.unmixing;
  cfg.topolabel = comp.topolabel;
  ica           = ft_componentanalysis(cfg, data);
  clear data
  
  % experimental conditions
  if isequal(tcgver, 1) % known vs. novel
    idx1            = find(ica.trialinfo(:,2)==1); % 1  = known data
    idx2            = find(ica.trialinfo(:,2)>1);  % 2-5 = novel data
  elseif isequal(tcgver, 2) % child vs. adult
    idx1            = find(ica.trialinfo(:,3)==1); % 1 = child addressee
    idx2            = find(ica.trialinfo(:,3)==2); % 2 = adult addressee
  end
  
  % spectral analysis
  cfg           = [];
  cfg.method    = 'mtmfft';
  cfg.taper     = 'hanning';
  cfg.foilim    = [1 300];
  cfg.trials    = idx1;
  freq1         = ft_freqanalysis(cfg, ica);
  cfg.trials    = idx2;
  freq2         = ft_freqanalysis(cfg, ica);
  
  % ERP analysis
  cfg           = [];
  cfg.vartrllength = 2;
  cfg.trials    = idx1;
  tlock1        = ft_timelockanalysis(cfg, ica);
  cfg.trials    = idx2;
  tlock2        = ft_timelockanalysis(cfg, ica);
  cfg           = [];
  cfg.baseline  = [-.5 0];
  tlock1        = ft_timelockbaseline(cfg, tlock1);
  tlock2        = ft_timelockbaseline(cfg, tlock2);
  
  % visualize topographies, weights, powerspectra, and ERPs
  ecogidx = match_str(comp.topolabel, ft_channelselection([subj(tcgver).ecog{:}], comp.topolabel));
  seegidx = match_str(comp.topolabel, ft_channelselection([subj(tcgver).seeg{:}], comp.topolabel));
  comp.topo       = abs(comp.topo); % weight signs are ambiguous
  u               = 5; % plots per figure
  for s = 1:ceil(numel(comp.label)/u) % for each figure set
    if s*u < numel(comp.label)
      sel = s*u-u+1:s*u;
    else
      sel = s*u-u+1:numel(comp.label);
      u = sel(end)-sel(1)+1;
    end
    figure; hold on
    for p = 1:u
      
      % topography
      subplot(u, 4, p*4-3)
      cfg             = [];
      cfg.component   = sel(p);
      cfg.layout      = lay;
      cfg.showoutline = 'yes';
      cfg.zlim        = 'zeromax';
      ft_topoplotIC(cfg, comp);
      
      % channel weights
      subplot(u, 4, p*4-2) % larger p-value ~ spatially more broad
      if ~isempty(ecogidx)
        hold on; plot(ecogidx, comp.topo(ecogidx,sel(p)), 'k');
      end
      if ~isempty(seegidx)
        hold on; plot(seegidx, comp.topo(seegidx,sel(p)), 'g');
      end
      title(['p = ' num2str(round(comp.prob(sel(p)),2))])
      xlim([1 numel(comp.label)])
      
      % powerspectra
      subplot(u, 4, p*4-1)
      loglog(freq1.freq, mean(freq1.powspctrm,1), 'k', 'linewidth', .5)
      hold on; loglog(freq1.freq, freq1.powspctrm(sel(p),:), 'b', 'linewidth', 1)
      hold on; loglog(freq2.freq, freq2.powspctrm(sel(p),:), 'r', 'linewidth', 1)
      xlim([1 300])
      
      % ERPs
      subplot(u, 4, p*4)
      plot(tlock1.time, tlock1.avg(sel(p),:), 'b', 'linewidth', 1)
      hold on; plot(tlock2.time, tlock2.avg(sel(p),:), 'r', 'linewidth', 1)
      xlim([-.5 1])
    end
    
    print([subj(tcgver).datadir 'ica' filesep subjectm(9:end) '_comp_' num2str(tcgver) '_' num2str(sel(1)) '-' num2str(sel(end)) '.pdf'], '-dpdf');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % increase for visibility
  end
  
  % databrowser
  cfg           = [];
  cfg.viewmode  = 'component'; 
  cfg.layout    = lay; 
  ft_databrowser(cfg, ica);
  
  % write down the bad components
  rejcomp = input('specify components for rejection (e.g. [1 2 6]): ');
  close all
  
  % enter visually detected artifacts in subject m-file
  fid = fopen([subjectm '.m'],'At');
  fprintf(fid,'\n%s\n',['%% Entered @ ' datestr(now)]);
  fprintf(fid,'%s\n', ['subj(' num2str(tcgver) ').badcomp = [' num2str(rejcomp) '];']);
  fclose all;
  
end
