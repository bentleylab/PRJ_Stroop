function func_envelopes_oneFband_resamp1K(PRCSDATADIR, DATASET, FREQBAND)%, PERIOD, TIMING)
% For all frequencies, windows, and elecs, calculate envelope means across all time

%% Set up variables and directories
% EEGLab won't call Octave functions in MATLAB, so load specific folders or
%       genpath and then rmpath the Octave functions
addpath('/home/knight/duration_PLV/Scripts/');

PrcsDataDir = PRCSDATADIR;
uScorePos   = strfind(DATASET,'_');% periodPos = strfind(DATASET,'.');
SBJ         = DATASET(1:uScorePos-1);
task        = DATASET(uScorePos+1:end);

%% Set Frequency and Timing parameters
% Load necessary variables
load(strcat(PrcsDataDir,'valid_elecs_allSBJandTasks.mat'));
validE = eval(['valid_elecs.' SBJ '.' task ';']);
SBJDataDir=strcat(PrcsDataDir,SBJ,'/',task,'/');
load(strcat(SBJDataDir,'gdat_notch.mat'));
load(strcat(SBJDataDir,'subj_globals.mat'),'srate');

% fBandListName=strcat(SBJDataDir,'fBand_Detection/',SBJ,'_',task,'_fBand_list.mat');
% load(fBandListName);
% FREQBAND = 'HG';

% Prepare Directories of Results/Output
envelopeDir = strcat(SBJDataDir,'envelopes/');
if ~isequal(exist(envelopeDir, 'dir'),7)
    mkdir(envelopeDir);
end

%% PLV Calculations
% Initialize matrix for all fbands and elecs
if srate~=1000
    [p,q] = rat(1000/srate);
    ts_length = length(resample(gdat(1,:), p, q));
else
    ts_length = length(gdat(1,:));
end
envelopes = NaN([length(validE) ts_length]);

for elec_idx = 1:length(validE)
    curr_elec = validE(elec_idx);
    disp(strcat('Calculating envelopes for elec=',num2str(curr_elec),...
        ' (',num2str(elec_idx),'/',num2str(length(validE)),')'));
    
    % Resample data to 1000 Hz
    if srate~=1000
        [p,q] = rat(1000/srate);
        elec_data = resample(gdat(curr_elec,:), p, q); % Careful this isn't using EEGlab firls
    else
        elec_data = gdat(curr_elec,:);
    end
    
    % Filter elec data
    bp_limits = func_return_bp_limits(SBJ,task,FREQBAND);
    disp(strcat(FREQBAND,'[lo,hi] = ',num2str(bp_limits(1)),'-',num2str(bp_limits(2))));
    elec_bandpass = func_EEGlab_bandpass(elec_data, srate, bp_limits(1), bp_limits(2));
    envelopes(elec_idx,:) = abs(hilbert(elec_bandpass));
end

%% Save envelope mean trial matrices
env_filename = strcat(envelopeDir,'envelopes_',FREQBAND,'_1kHz.mat');
save(env_filename,'envelopes');

end
