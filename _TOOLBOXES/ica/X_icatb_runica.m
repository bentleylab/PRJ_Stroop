function [weights,sphere,activations,bias,signs,lrates,y] = icatb_runica(data,varargin)
% runica() - Perform Independent Component Analysis (ICA) decomposition
%                of psychophysiological data using the infomax ICA algorithm of 
%                Bell & Sejnowski (1995) with the natural gradient feature 
%                of Amari, Cichocki & Yang, the extended-ICA algorithm 
%                of Lee, Girolami & Sejnowski, PCA dimension reduction,
%                and/or specgram() preprocessing (M. Zibulevsky).
% Usage:
%        simply >> [weights,sphere] = runica(data);
%         or     
%          else >> [weights,sphere,activations,bias,signs,lrates] ...
%                                            = runica(data,'Key1',Value1',...);
% Input_Variable:
%
% data          = input data (chans,frames*epochs). 
%                    Note: If data consists of multiple discontinuous epochs, 
%                    each epoch should be separately baseline-zero'd using:
%                        >> data = rmbase(data,frames,basevector);
%
% Optional_Keywords         Keyword_Values                        Default_Values
%
% 'ncomps'     = [N] number of ICA components to compute (default -> chans)
%                    using rectangular ICA decomposition
% 'pca'         = [N] decompose a principal component      (default -> 0=off)
%                    subspace of the data. Value is the number of PCs to retain.
% 'sphering'  = ['on'/'off'] flag sphering of data        (default -> 'on')
% 'weights'    = [W] initial weight matrix                    (default -> eye())
%                                     (Note: if 'sphering' 'off', default -> spher())
% 'lrate'      = [rate] initial ICA learning rate (<< 1) (default -> heuristic)
% 'block'      = [N] ICA block size (<< datalength)        (default -> heuristic)
% 'anneal'     = annealing constant (0,1] (defaults -> 0.90, or 0.98, extended)
%                                 controls speed of convergence
% 'annealdeg' = [N] degrees weight change for annealing (default -> 70)
% 'stop'        = [f] stop training when weight-change < this (default -> 1e-6)
% 'maxsteps'  = [N] max number of ICA training steps     (default -> 512)
% 'bias'        = ['on'/'off'] perform bias adjustment     (default -> 'on')
% 'momentum'  = [0<f<1] training momentum                    (default -> 0)
% 'extended'  = [N] perform tanh() "extended-ICA" with sign estimation 
%                    every N training blocks. If N < 0, fix number of sub-Gaussian
%                    components to -N [faster than N>0]        (default|0 -> off)
% 'specgram'  = [srate loHz hiHz Hzinc frames] decompose a complex time/frequency
%                    transform of the data            (defaults [srate 0 srate/2 1 0])
% 'posact'     = make all component activations net-positive(default 'on')
% 'verbose'    = give ascii messages ('on'/'off')          (default -> 'on')
% 'gpu'        = gpu support with MATLAB Parallel Computing Toolbox and Jacket
%                    GPU Engine('pct'/'jkt'/'off')         (default -> 'off')
%
% Output_Variables [RO = output in reverse order of projected mean variance 
%                                unless starting weight matrix passed ('weights' above)]
%
% weights      = ICA weight matrix (comps,chans)      [RO]
% sphere        = data sphering matrix (chans,chans) = spher(data)
%                    Note: unmixing_matrix = weights*sphere {sphering off -> eye(chans)}
% activations = activation time courses of the output components (ncomps,frames*epochs)
% bias          = vector of final (ncomps) online bias [RO]     (default = zeros())
% signs         = extended-ICA signs for components     [RO]     (default = ones())
%                         [-1 = sub-Gaussian; 1 = super-Gaussian]
% lrates        = vector of learning rates used at each training step
%

% Toolbox Citation:
%
% Makeig, Scott et al. "ICA Toolbox for Psychophysiological Research (version 3.4)". 
% WWW Site, Computational Neurobiology Laboratory, The Salk Institute for Biological 
% Studies <www.cnl.salk.edu/~ica.html>, 1999. [World Wide Web Publication]. 
%
% 1st Publication:
%
% Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J.,
% "Independent component analysis of electroencephalographic data," 
% In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural 
% Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).
%
% For more information:
% http://www.cnl.salk.edu/~scott/icafaq.html - FAQ on ICA/EEG
% http://www.cnl.salk.edu/~scott/icabib.html - mss. on ICA & biosignals
% http://www.cnl.salk.edu/~tony/ica.html - math. mss. on ICA, with kernal code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit history %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  runica()  - by Scott Makeig with contributions from Tony Bell, Te-Won Lee 
%                  Tzyy-Ping Jung, Sigurd Enghoff, Michael Zibulevsky et al.
%                                     CNL / Salk Institute 1996-99
%  04-30-96 built from icatest.m and ~jung/.../wtwpwica.m -sm
%  07-28-97 new runica(), adds bias (default on), momentum (default off), 
%              extended-ICA (Lee & Sejnowski, 1997), cumulative angledelta 
%              (until lrate drops), keywords, signcount for speeding extended-ICA
%  10-07-97 put acos() outside verbose loop; verbose 'off' wasn't stopping -sm
%  11-11-97 adjusted help msg -sm
%  11-30-97 return eye(chans) if sphering 'off' or 'none' (undocumented option) -sm
%  02-27-98 use pinv() instead of inv() to rank order comps if ncomps < chans -sm
%  04-28-98 added 'posact' and 'pca' flags  -sm
%  07-16-98 reduced length of randperm() for kurtosis subset calc. -se & sm
%  07-19-98 fixed typo in weights def. above -tl & sm
%  12-21-99 added 'specgram' option suggested by Michael Zibulevsky, UNM -sm
%  12-22-99 fixed rand() sizing inefficiency on suggestion of Mike Spratling, UK -sm
%  04-27-12 added GPU computing support, NIH/NIDCD -Yisheng Xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    help runica  
    return
end

[chans frames] = size(data); % determine the data size
urchans = chans;  % remember original data channels 
datalength = frames;
if chans<2 || frames<chans
  fprintf('\nrunica() - data size (%d,%d) too small.\n\n', chans,frames);
  return
end

%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
%
MAX_WEIGHT          = 1e8;          % guess that weights larger than this have blown up
DEFAULT_STOP        = 0.000001;     % stop training if weight changes below this
DEFAULT_ANNEALDEG   = 60;           % when angle change reaches this value,
DEFAULT_ANNEALSTEP  = 0.90;         %     anneal by multiplying lrate by this
DEFAULT_EXTANNEAL   = 0.98;         %     or this if extended-ICA
DEFAULT_MAXSTEPS    = 512;          % ]top training after this many steps 
DEFAULT_MOMENTUM    = 0.0;          % default momentum weight

DEFAULT_BLOWUP      = 1000000000.0; % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC  = 0.8;          % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC = 0.9;          % if weights blowup, restart with lrate
% lower by this factor
MIN_LRATE           = 0.000001;     % if weight blowups make lrate < this, quit
MAX_LRATE           = 0.1;          % guard against uselessly high learning rate
DEFAULT_LRATE       = 0.015/log(chans); 
% heuristic default - may need adjustment
%    for large or tiny data sets!
DEFAULT_BLOCK       = floor(sqrt(frames/3));    % heuristic default
% - may need adjustment!
% Extended-ICA option:
DEFAULT_EXTENDED    = 0;            % default off
DEFAULT_EXTBLOCKS   = 1;            % number of blocks per kurtosis calculation
DEFAULT_NSUB        = 1;            % initial default number of assumed sub-Gaussians
% for extended-ICA
DEFAULT_EXTMOMENTUM = 0.5;          % momentum term for computing extended-ICA kurtosis
MAX_KURTSIZE        = 6000;         % max points to use in kurtosis calculation
MIN_KURTSIZE        = 2000;         % minimum good kurtosis size (flag warning)
SIGNCOUNT_THRESHOLD = 25;           % raise extblocks when sign vector unchanged
% after this many steps
SIGNCOUNT_STEP      = 2;            % extblocks increment factor 

DEFAULT_SPHEREFLAG  = 'on';         % use the sphere matrix as the default
%    starting weight matrix
DEFAULT_PCAFLAG     = 'off';        % don't use PCA reduction
DEFAULT_POSACTFLAG  = 'on';         % use posact()
DEFAULT_VERBOSE     = 1;            % write ascii info to calling screen
DEFAULT_BIASFLAG    = 1;            % default to using bias in the ICA update rule
DEFAULT_GPU         = 3;            % default to using CPU
GPU_ENGINES         = {'pct','jkt','off'};      % GPU computing engines
GPU_DESCRIPTION     = {'with MATLAB PCT/GPU ','with Jacket GPU Engine ',''};
GPU_CAST_FUN        = {'gpuArray(','(g','('};   % CPU-GPU cast functions
CPU_CAST_FUN        = {'gather','double',''};   % GPU-CPU cast functions
%                                            
%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargout < 2, 
    fprintf('runica() - needs at least two output arguments.\n');
    return
end
epochs = 1;							% do not care how many epochs in data

pcaflag     = DEFAULT_PCAFLAG;
sphering    = DEFAULT_SPHEREFLAG;   % default flags
posactflag  = DEFAULT_POSACTFLAG;
verbose     = DEFAULT_VERBOSE;

block       = DEFAULT_BLOCK;        % heuristic default - may need adjustment!
lrate       = DEFAULT_LRATE;
annealdeg   = DEFAULT_ANNEALDEG;
annealstep  = 0;                    % defaults declared below
nochange    = DEFAULT_STOP;
momentum    = DEFAULT_MOMENTUM;
maxsteps    = DEFAULT_MAXSTEPS;

weights     = 0;                    % defaults defined below
ncomps      = chans;
biasflag    = DEFAULT_BIASFLAG;

extended    = DEFAULT_EXTENDED;
extblocks   = DEFAULT_EXTBLOCKS;
kurtsize    = MAX_KURTSIZE;
signsbias   = 0.02;                 % bias towards super-Gaussian components
extmomentum = DEFAULT_EXTMOMENTUM;  % exp. average the kurtosis estimates
nsub        = DEFAULT_NSUB;
wts_blowup  = 0;                    % flag =1 when weights too large
wts_passed  = 0;                    % flag weights passed as argument
gpu         = DEFAULT_GPU;
%
%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%
if (nargin > 1 && rem(nargin,2) == 0)
    fprintf('runica(): Even number of input arguments???')
    return
end
for i = 1:2:nargin-2    % for each Keyword
    Keyword = varargin{i};
    Value = varargin{i+1};
    if ~ischar(Keyword)
        fprintf('runica(): keywords must be strings')
        return
    end
    
    switch lower(Keyword)
      case {'weights','weight'}
        if ischar(Value)
            fprintf(...
                'runica(): weights value must be a weight matrix or sphere')
            return
        else
            weights = Value;
            wts_passed =1;
        end
      case 'ncomps'
        if ischar(Value)
            fprintf('runica(): ncomps value must be an integer')
            return
        end
        if ncomps < urchans && ncomps ~= Value
            fprintf('runica(): Use either PCA or ICA dimension reduction');
            return
        end
        ncomps = Value;
        if ~ncomps,
            ncomps = chans;
        end
      case 'pca'
        if ncomps < urchans && ncomps ~= Value
            fprintf('runica(): Use either PCA or ICA dimension reduction');
            return
        end
        if ischar(Value)
            fprintf(...
                'runica(): pca value should be the number of principal components to retain')
            return
        end
        pcaflag = 'on';
        ncomps = Value;
        if ncomps >= chans || ncomps < 1,
            fprintf('runica(): pca value must be in range [1,%d]\n',chans-1)
            return
        end
        chans = ncomps;
      case 'posact'
        if ~ischar(Value)
            fprintf('runica(): posact value must be on or off')
            return
        else 
            Value = lower(Value);
            if ~strcmp(Value,'on') && ~strcmp(Value,'off'),
                fprintf('runica(): posact value must be on or off')
                return
            end
            posactflag = Value;
        end
      case 'lrate'
        if ischar(Value)
            fprintf('runica(): lrate value must be a number')
            return
        end
        lrate = Value;
        if lrate > MAX_LRATE || lrate < 0,
            fprintf('runica(): lrate value is out of bounds'); 
            return
        end
        if ~lrate,
            lrate = DEFAULT_LRATE;
        end
      case {'block','blocksize'}
        if ischar(Value)
            fprintf('runica(): block size value must be a number')
            return
        end
        block = Value;
        if ~block,
            block = DEFAULT_BLOCK; 
        end
      case {'stop','nochange','stopping'}
        if ischar(Value)
            fprintf('runica(): stop wchange value must be a number')
            return
        end
        nochange = Value;
      case {'maxsteps','steps'}
        if ischar(Value)
            fprintf('runica(): maxsteps value must be an integer')
            return
        end
        maxsteps = Value;
        if ~maxsteps,
            maxsteps    = DEFAULT_MAXSTEPS;
        end
        if maxsteps < 0
            fprintf('runica(): maxsteps value must be a positive integer')
            return
        end
      case {'anneal','annealstep'}
        if ischar(Value)
            fprintf('runica(): anneal step constant must be a number (0,1)')
            return
        end
        annealstep = Value;
        if annealstep <=0 || annealstep > 1,
            fprintf('runica(): anneal step value must be (0,1]')
            return
        end
      case {'annealdeg','degrees'}
        if ischar(Value)
            fprintf('runica(): annealdeg value must be a number')
            return
        end
        annealdeg = Value;
        if ~annealdeg,
            annealdeg = DEFAULT_ANNEALDEG;
        elseif annealdeg > 180 || annealdeg < 0
            fprintf('runica(): annealdeg value is out of bounds [0,180]')
            return
            
        end
      case 'momentum'
        if ischar(Value)
            fprintf('runica(): momentum value must be a number')
            return
        end
        momentum = Value;
        if momentum > 1.0 || momentum < 0
            fprintf('runica(): momentum value is out of bounds [0,1]')
            return
        end
      case {'sphering','sphereing','sphere'}
        if ~ischar(Value)
            fprintf('runica(): sphering value must be on, off, or none')
            return
        else 
            Value = lower(Value);
            if ~strcmp(Value,'on') && ~strcmp(Value,'off') && ~strcmp(Value,'none'),
                fprintf('runica(): sphering value must be on or off')
                return
            end
            sphering = Value;
        end
      case 'bias'
        if ~ischar(Value)
            fprintf('runica(): bias value must be on or off')
            return
        else 
            Value = lower(Value);
            if strcmp(Value,'on') 
                biasflag = 1;
            elseif strcmp(Value,'off'),
                biasflag = 0;
            else
                fprintf('runica(): bias value must be on or off')
                return
            end
        end
      case {'specgram','spec'}
        if ischar(Value)
            fprintf('runica(): specgram argument must be a vector')
            return
        end
        srate = Value(1);
        if (srate < 0)
            fprintf('runica(): specgram srate must be >=0')
            return
        end
        if length(Value)>1
            loHz = Value(2);
            if (loHz < 0 || loHz > srate/2)
                fprintf('runica(): specgram loHz must be >=0 and <= srate/2')
                return
            end
        else
            loHz = 0;       % default
        end
        if length(Value)>2
            hiHz = Value(3);
            if (hiHz < loHz || hiHz > srate/2)
                fprintf('runica(): specgram hiHz must be >=loHz and <= srate/2')
                return
            end
        else
            hiHz = srate/2; % default
        end
        if length(Value)>3
            Hzinc = Value(4);
            if (Hzinc<=0 || Hzinc>hiHz-loHz)
                fprintf('runica(): specgram Hzinc must be >0 and <= hiHz-loHz')
                return
            end
        else
            Hzinc = 1;      % default
        end
        if length(Value)>4
            Hzframes = Value(5);
            if (Hzframes<0 || Hzframes > size(data,2))
                fprintf('runica(): specgram frames must be >=0 and <= data length')
                return
            end
        else
            Hzframes = size(data,2);    % default
        end
      case {'extended','extend'}
        if ischar(Value)
            fprintf('runica(): extended value must be an integer (+/-)')
            return
        else
            extended = 1;           % turn on extended-ICA
            extblocks = fix(Value); % number of blocks per kurt() compute
            if extblocks < 0
                nsub = -1*fix(extblocks);   % fix this many sub-Gauss comps
            elseif ~extblocks,
                extended = 0;               % turn extended-ICA off
            elseif kurtsize>frames,         % length of kurtosis calculation
                kurtsize = frames;
                if kurtsize < MIN_KURTSIZE
                    fprintf(...
                        'runica() warning: kurtosis values inexact for << %d points.\n',...
                        MIN_KURTSIZE);
                end
            end
        end
      case 'verbose'
        if ~ischar(Value)
            fprintf('runica(): verbose flag value must be on or off')
            return
        else
            Value = lower(Value);
            if strcmp(Value,'on'),
                verbose = 1; 
            elseif strcmp(Value,'off'),
                verbose = 0; 
            else
                fprintf('runica(): verbose flag value must be on or off')
                return
            end
        end
      case 'gpu'
        if ~ischar(Value)
            fprintf('runica(): gpu flag value must be a string')
            return
        else
            gpu = find(strcmpi(Value,GPU_ENGINES));
            if gpu == 0,
                fprintf('runica(): unkown gpu flag value')
                return
            end
        end
      otherwise
        fprintf('runica(): unknown flag')
        return
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize weights, etc. %%%%%%%%%%%%%%%%%%%%%%%%
%
if ~annealstep,
    if ~extended,
        annealstep = DEFAULT_ANNEALSTEP;        % defaults defined above
    else
        annealstep = DEFAULT_EXTANNEAL;         % defaults defined above
    end
end % else use annealstep from commandline

if ~annealdeg, 
    annealdeg  = DEFAULT_ANNEALDEG - momentum*90;   % heuristic
    if annealdeg < 0,
        annealdeg = 0;
    end
end
if ncomps >  chans || ncomps < 1
    fprintf('runica(): number of components must be 1 to %d.\n',chans);
    return
end

if weights ~= 0,                                % initialize weights
    % starting weights are being passed to runica() from the commandline
    if verbose,
        fprintf('Using starting weight matrix named in argument list ...\n')
    end
    if  chans>ncomps && weights ~=0,
        [r,c]=size(weights);
        if r~=ncomps || c~=chans,
            fprintf(...
                'runica(): weight matrix must have %d rows, %d columns.\n', ...
                chans,ncomps);
            return;
        end
    end
end;    
%
%%%%%%%%%%%%%%%%%%%%% Check keyword values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if frames<chans,
    fprintf('runica(): data length %d < data channels %f!\n',frames,chans)
    return
elseif block < 2,
    fprintf('runica(): block size %d too small!\n',block)
    return
elseif block > frames, 
    fprintf('runica(): block size exceeds data length!\n');
    return
elseif floor(epochs) ~= epochs,
    fprintf('runica(): data length is not a multiple of the epoch length!\n');
    return
elseif nsub > ncomps
    fprintf('runica(): there can be at most %d sub-Gaussian components!\n',ncomps);
    return
end;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the data %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf( ...
        '\nInput data size [%d,%d] = %d channels, %d frames.\n', ...
        chans,frames,chans,frames);
    if strcmp(pcaflag,'on')
        fprintf('After PCA dimension reduction,\n  finding ');
    else
        fprintf('Finding ');
    end
    if ~extended
        fprintf('%d ICA components using logistic ICA.\n',ncomps);
    else    % if extended
        fprintf('%d ICA components using extended ICA.\n',ncomps);
        if extblocks > 0
            fprintf(...
                'Kurtosis will be calculated initially every %d blocks using %d data points.\n',...
                extblocks,      kurtsize);
        else
            fprintf(...
                'Kurtosis will not be calculated. Exactly %d sub-Gaussian components assumed.\n',...
                nsub);
        end
    end
    fprintf('Initial learning rate will be %g, block size %d.\n',lrate,block);
    if momentum>0,
        fprintf('Momentum will be %g.\n',momentum);
    end
    fprintf( ...
        'Learning rate will be multiplied by %g whenever angledelta >= %g deg.\n', ...
        annealstep,annealdeg);
    fprintf('Training will end when wchange < %g or after %d steps.\n', ...
        nochange,maxsteps);
    if biasflag,
        fprintf('Online bias adjustment will be used.\n');
    else
        fprintf('Online bias adjustment will not be used.\n');
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Remove overall row means %%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf('Removing mean of each channel ...\n');
end
 disp('Not removing mean of each channel!!!');
%data = data - mean(data')'*ones(1,frames);        % subtract row means

if verbose,
    fprintf('Final training data range: %g to %g\n', ...
        min(min(data)),max(max(data)));
end
%
%%%%%%%%%%%%%%%%%%% Perform PCA reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Reducing the data to %d principal dimensions...\n',ncomps);
    [eigenvectors,eigenvalues,data] = icatb_pcsquash(data,ncomps);
    % make data its projection onto the ncomps-dim principal subspace
end
%
%%%%%%%%%%%%%%%%%%% Perform specgram transformation %%%%%%%%%%%%%%%%%%%%%%%
%
if exist('srate','var') 
    % [P F T] = SPECGRAM(A,NFFT,Fs,WINDOW,NOVERLAP)
    Hzwinlen =  fix(srate/Hzinc);
    Hzfftlen = 2^(ceil(log(Hzwinlen)/log(2)));
    if (Hzwinlen>Hzframes)
        Hzwinlen = Hzframes;
    end
    Hzoverlap = 0;
    if (Hzwinlen>Hzframes)
        Hzoverlap = Hzframes - Hzwinlen;
    end
    % get freqs and times
    [tmp,freqs,tms] = specgram(data(1,:),Hzfftlen,srate,Hzwinlen,Hzoverlap);
    fs = find(freqs>=loHz & freqs <= hiHz);
    % fprintf('    size(fs) = %d,%d\n',size(fs,1),size(fs,2));
    % fprintf('    size(tmp) = %d,%d\n',size(tmp,1),size(tmp,2));
    speclength = length(fs)*size(tmp,2);
    datalength = speclength*2;
    specdata = zeros(chans,datalength);
    tmp = reshape(tmp(fs,:),1,speclength);
    specdata(1,:) = [real(tmp) imag(tmp)];
    for ch=2:chans
        tmp = specgram(data(ch,:),Hzwinlen,srate,Hzwinlen,Hzoverlap);
        tmp = reshape((tmp(fs,:)),1,speclength);
        specdata(ch,:) = [real(tmp) imag(tmp)]; % channels are rows
    end
    fprintf('Converted data to %d channels by %d=2*%dx%d points spectrogram data.\n',chans,2*length(fs)*length(tms),length(fs),length(tms));
    fprintf('    Low Hz %g, high Hz %g, Hz inc %g, window length %d\n',freqs(fs(1)),freqs(fs(end)),freqs(fs(2))-freqs(fs(1)),Hzwinlen);
    data = specdata;
end
%
%%%%%%%%%%%%%%%%%%% Perform sphering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if strcmp(sphering,'on'), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verbose,
        fprintf('Computing the sphering matrix...\n');
    end
    sphere = 2.0*inv(sqrtm(cov(data'))); %#ok<MINV> % find the "sphering" matrix = spher()
    if ~weights,
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
        end
        weights = eye(ncomps,chans); %#ok<*NASGU> % begin with the identity matrix
    else % weights given on commandline
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
        end
    end
    if verbose,
        fprintf('Sphering the data ...\n');
    end
    data = sphere*data;        % actually decorrelate the electrode signals
    
elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~weights
        if verbose,
            fprintf('Using the sphering matrix as the starting weight matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere = 2.0*inv(sqrtm(cov(data')));    %#ok<MINV> % find the "sphering" matrix = spher()
        weights = eye(ncomps,chans)*sphere;     % begin with the identity matrix
        sphere = eye(chans);                    % return the identity matrix
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        sphere = eye(chans);                    % return the identity matrix
    end
elseif strcmp(sphering,'none')
    sphere = eye(chans);                        % return the identity matrix
    if ~weights
        if verbose,
            fprintf('Starting weights are the identity matrix ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        weights = eye(ncomps,chans); % begin with the identity matrix
    else % weights ~= 0
        if verbose,
            fprintf('Using starting weights named on commandline ...\n');
            fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
    end
    if verbose,
        fprintf('Returned variable "sphere" will be the identity matrix.\n');
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
%
signs = ones(1,ncomps);         % initialize signs to nsub -1, rest +1
for k=1:nsub
    signs(k) = -1;
end
if extended && extblocks < 0 && verbose,
    fprintf('Fixed extended-ICA sign assignments:  ');
    for k=1:ncomps
        fprintf('%d ',signs(k));	
    end; fprintf('\n');
end
BI = eval(['block*' GPU_CAST_FUN{gpu} 'eye(ncomps,ncomps,''double''))']);
weights = eval([GPU_CAST_FUN{gpu} 'double(weights))']);
prevwtchange = eval([GPU_CAST_FUN{gpu} 'zeros(chans,ncomps,''double''))']);
onesrow = eval([GPU_CAST_FUN{gpu} 'ones(1,block,''double''))']);
bias = eval([GPU_CAST_FUN{gpu} 'zeros(ncomps,1,''double''))']);
signs = eval(['diag(' GPU_CAST_FUN{gpu} 'double(signs)))']);
oldsigns = eval([GPU_CAST_FUN{gpu} 'zeros(size(signs),''double''))']);
old_kk = eval([GPU_CAST_FUN{gpu} 'zeros(1,ncomps,''double''))']);
lastt = fix((datalength/block-1)*block+1);
changes = [];
degconst = 180./pi;
startweights = weights;
prevweights = startweights;
oldweights = startweights;
lrates = zeros(1,maxsteps);
signcount = 0;                  % counter for same-signs
signcounts = [];
urextblocks = extblocks;        % original value, for resets
%
%%%%%%%% ICA training loop using the logistic sigmoid %%%%%%%%%%%%%%%%%%%
%
if verbose,
    fprintf('Beginning ICA training %s...',GPU_DESCRIPTION{gpu});
    if extended,
        fprintf(' first training step may be slow ...\n');
    else
        fprintf('\n');
    end
end
step=0;
laststep=0; 
blockno = 1;  % running block counter for kurtosis interrupts

while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    permute=randperm(datalength); % shuffle data order at each step
    
    for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
        if biasflag                                                                    
            u=weights*data(:,permute(t:t+block-1)) + bias*onesrow;        
        else                                                                                 
            u=weights*data(:,permute(t:t+block-1));                             
        end                                                                                  
        if ~extended
            %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
            y=1./(1+exp(-u));                                                %
            weights=weights+lrate*(BI+(1-2*y)*u')*weights;                   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else % extended-ICA
            %%%%%%%%%%%%%%%%%%% Extended-ICA weight update %%%%%%%%%%%%%%%%%%%
            y=tanh(u);                                                       %
            weights = weights + lrate*(BI-signs*y*u'-u*u')*weights;          %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if biasflag 
            if ~extended
                %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
                bias = bias + lrate*sum(1-2*y,2);   % for logistic nonlin.       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % extended
                %%%%%%%%%%%%%%%%%%% Extended-ICA bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                bias = bias + lrate*sum(-2*y,2);    % for tanh() nonlin.         %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end                                                
        end
        
        if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            weights = weights + momentum*prevwtchange;                     
            prevwtchange = weights-prevweights;                             
            prevweights = weights;                                             
        end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if eval(['max(max(' CPU_CAST_FUN{gpu} '(abs(weights))))']) > MAX_WEIGHT
            wts_blowup = 1;
            change = nochange;
        end
        if extended && ~wts_blowup
            %
            %%%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
            %
            if extblocks > 0 && rem(blockno,extblocks) == 0, 
                % recompute signs vector using kurtosis
                if kurtsize < frames 
                    %rp = ceil(rand(kurtsize)*datalength);     % pick random subset
                    rp = randperm(datalength);
                    % of data. Compute
                    partact=weights*data(:,rp(1:kurtsize)); % partial activation
                else                                                     % for small data sets,
                    partact=weights*data;                         % use whole data
                end
                m2=mean(partact'.^2).^2; 
                m4= mean(partact'.^4);
                kk= (m4./m2)-3.0;                                    % kurtosis estimates
                if extmomentum
                    kk = extmomentum*old_kk + (1.0-extmomentum)*kk; % use momentum
                    old_kk = kk;
                end
                signs = diag(sign(kk+signsbias));  % pick component signs
                if eval([CPU_CAST_FUN{gpu} '(signs) == ' CPU_CAST_FUN{gpu} '(oldsigns)']),
                    signcount = signcount+1;
                else
                    signcount = 0;
                end
                oldsigns = signs;
                signcounts = [signcounts signcount]; %#ok<AGROW>
                if signcount >= SIGNCOUNT_THRESHOLD,
                    extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
                    signcount = 0;                                      % less frequent if sign
                end                                                      % is not changing
            end % extblocks > 0 & . . .
        end % if extended %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        blockno = blockno + 1;
        if wts_blowup
            break
        end
    end % training block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~wts_blowup
        oldwtchange = weights-oldweights;
        step=step+1; 
        %
        %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
        %
        lrates(1,step) = lrate;
        angledelta = 0;
        delta = reshape(oldwtchange,1,chans*ncomps);
        change = eval([CPU_CAST_FUN{gpu} '(delta*delta'')']); 
    end
    %
    %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
    %
    if wts_blowup || isnan(change) || isinf(change),  % if weights blow up,
        fprintf('');
        bias = eval([GPU_CAST_FUN{gpu} 'zeros(ncomps,1,''double''))']);
        prevwtchange = eval([GPU_CAST_FUN{gpu} 'zeros(chans,ncomps,''double''))']);
        delta = eval([GPU_CAST_FUN{gpu} 'zeros(1,chans*ncomps,''double''))']);
        olddelta = delta;
        step = 0;                                  % start again
        wts_blowup = 0;                          % re-initialize variables
        blockno = 1;
        lrate = lrate*DEFAULT_RESTART_FAC; % with lower learning rate
        weights = startweights;                % and original weight matrix
        oldweights = startweights;                
        prevweights = startweights;
        change = nochange;
        extblocks = urextblocks;
        lrates = zeros(1,maxsteps);
        if extended
            signs = ones(1,ncomps); % initialize signs to nsub -1, rest +1
            for k=1:nsub
                signs(k) = -1;
            end
            signs = eval(['diag(' GPU_CAST_FUN{gpu} 'double(signs)))']);
            oldsigns = eval([GPU_CAST_FUN{gpu} 'zeros(size(signs),''double''))']);
        end
        if lrate > MIN_LRATE
            r = rank(data);
            if r<ncomps
                fprintf('Data has rank %d. Cannot compute %d components.\n',...
                    r,ncomps);
                return
            else
                fprintf(...
                    'Lowering learning rate to %g and starting again.\n',lrate);
            end
        else
            fprintf( ...
                'runica(): QUITTING - weight matrix may not be invertible!\n');
            return;
        end
    else % if weights in bounds 
        %
        %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
        %
        if step > 2 
            angledelta = eval(['acos(' CPU_CAST_FUN{gpu} '(delta*olddelta'')/sqrt(change*oldchange))']);
        end
        if verbose,
            if step > 2, 
                if ~extended,
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
                        step,lrate,change,degconst*angledelta);
                else
                    fprintf(...
                        'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg, %d subgauss\n',...
                        step,lrate,change,degconst*angledelta,(ncomps-sum(diag(signs)))/2);
                end
            elseif ~extended
                fprintf(...
                    'step %d - lrate %5f, wchange %7.6f\n',step,lrate,change);
            else
                fprintf(...
                    'step %d - lrate %5f, wchange %7.6f, %d subgauss\n',...
                    step,lrate,change,(ncomps-sum(diag(signs)))/2);
            end % step > 2
        end; % if verbose
        %
        %%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        changes = [changes change]; %#ok<AGROW>
        oldweights = weights;
        %
        %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if degconst*angledelta > annealdeg,  
            lrate = lrate*annealstep;             % anneal learning rate
            olddelta    = delta;                     % accumulate angledelta until
            oldchange  = change;                    %  annealdeg is reached
        elseif step == 1                            % on first step only
            olddelta    = delta;                     % initialize 
            oldchange  = change;                    
        end
        %
        %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if step > 2 && change < nochange,        % apply stopping rule
            laststep = step;                
            step = maxsteps;                        % stop when weights stabilize
        elseif change > DEFAULT_BLOWUP,        % if weights blow up,
            lrate = lrate*DEFAULT_BLOWUP_FAC;     % keep trying 
        end;                                            % with a smaller learning rate
    end; % end if weights in bounds
    
end; % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weights = eval([CPU_CAST_FUN{gpu} '(weights)']);
if ~laststep
    laststep = step;
end;
lrates = lrates(1,1:laststep);              % truncate lrate history vector
%
%%%%%%%%%%%%%% Orient components towards positive activation %%%%%%%%%%%
%
if strcmp(posactflag,'on')
    [activations, winvout, weights] = icatb_posact(data, weights);
    % changes signs of activations and weights to make activations
    % net rms-positive
else
    activations = weights*data;
end
%
%%%%%%%%%%%%%% If pcaflag, compose PCA and ICA matrices %%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    fprintf('Composing the eigenvector, weights, and sphere matrices\n');
    fprintf('  into a single rectangular weights matrix; sphere=eye(%d)\n'...
        ,chans);
    weights= weights*sphere*eigenvectors(:,1:ncomps)'; 
    sphere = eye(urchans);
end
%
%%%%%% Sort components in descending order of max projected variance %%%%
%
if verbose,
    fprintf(...
        'Sorting components in descending order of mean projected variance ...\n');
end
if wts_passed == 0
    %
    %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    meanvar  = zeros(ncomps,1);        % size of the projections
    if ncomps == urchans % if weights are square . . .
        winv = inv(weights*sphere);
    else
        fprintf('Using pseudo-inverse of weight matrix to rank order component projections.\n');
        winv = pinv(weights*sphere);
    end
    for s=1:ncomps
        if verbose,
            fprintf('%d ',s);            % construct single-component data matrix
        end
        % project to scalp, then add row means 
        compproj = winv(:,s)*activations(s,:);
        meanvar(s) = mean(sum(compproj.*compproj)/(size(compproj,1)-1));
        % compute mean variance 
    end                                                      % at all scalp channels
    if verbose,
        fprintf('\n');
    end
    %
    %%%%%%%%%%%%%% Sort components by mean variance %%%%%%%%%%%%%%%%%%%%%%%%
    %
    [sortvar, windex] = sort(meanvar);
    windex = windex(ncomps:-1:1); % order large to small 
%     meanvar = meanvar(windex);
    % 
    %%%%%%%%%%%%%%%%%%%%% Filter data using final weights %%%%%%%%%%%%%%%%%%
    %
    if nargout > 2, % if activations are to be returned
        if verbose,
            fprintf('Permuting the activation wave forms ...\n');
        end
        activations = activations(windex,:);
    else
        clear activations
    end
    weights = weights(windex,:);% reorder the weight matrix
    bias  = eval([CPU_CAST_FUN{gpu} '(bias(windex))']);		% reorder them
    signs = diag(signs);          % vectorize the signs matrix
    signs = eval([CPU_CAST_FUN{gpu} '(signs(windex))']);    % reorder them
else
    fprintf('Components not ordered by variance.\n');
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if nargout > 6
    u=weights*data + bias*ones(1,frames);        
    y = zeros(size(u));
    for c=1:chans
        for f=1:frames
            y(c,f) = 1/(1+exp(-u(c,f)));
        end
    end
end
