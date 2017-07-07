function func_envelopes_allOpts_v2(PRCSDATADIR, DATASET)%, FREQBAND, PERIOD, TIMING)
% For all frequencies, windows, and elecs, calculate envelope means across all time

%% Set up variables and directories
% EEGLab won't call Octave functions in MATLAB, so load specific folders or
%       genpath and then rmpath the Octave functions
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/EEGLab/eeglab12_0_2_6b/'));

prcs_data_dir = PRCSDATADIR;
uScorePos   = strfind(DATASET,'_');% periodPos = strfind(DATASET,'.');
SBJ         = DATASET(1:uScorePos-1);
task        = DATASET(uScorePos+1:end);
WL          = 500;          % Window length for PLV in ms

%% Set Frequency and Timing parameters
% Load necessary variables
load(strcat(prcs_data_dir,'valid_elecs_allSBJandTasks.mat'));
validE = eval(['valid_elecs.' SBJ '.' task ';']);
SBJDataDir=strcat(prcs_data_dir,SBJ,'/',task,'/');
load(strcat(SBJDataDir,'gdat_notch.mat'));
load(strcat(SBJDataDir,'subj_globals.mat'),'srate');
WLinDP = WL*srate/1000;

fBandListName=strcat(SBJDataDir,'fBand_Detection/',SBJ,'_',task,'_fBand_list.mat');
load(fBandListName);
FREQBAND = {FREQBAND{:} 'HG'};

% Prepare Directories of Results/Output
envelopeDir = strcat(SBJDataDir,'envelopes/');
if ~isequal(exist(envelopeDir, 'dir'),7)
    mkdir(envelopeDir);
end

%% PLV Calculations
% Initialize matrix for all fbands and elecs
master_envelopes = NaN([length(FREQBAND) length(validE) size(gdat,2)]);

for elec_idx = 1:length(validE)
    curr_elec = validE(elec_idx);
    disp(strcat('Calculating envelopes for elec=',num2str(curr_elec),...
        ' (',num2str(elec_idx),'/',num2str(length(validE)),')'));
    
    for freqband = 1:length(FREQBAND)
        elec_data = gdat(curr_elec,:);
        
        % Filter elec data
        bp_limits = func_return_bp_limits(SBJ,task,FREQBAND{freqband});
        disp(strcat(FREQBAND{freqband},'[lo,hi] = ',num2str(bp_limits(1)),'-',num2str(bp_limits(2))));
        elec_highpass = eegfilt(elec_data, srate, bp_limits(1), []);
        elec_bandpass = eegfilt(elec_highpass, srate, [], bp_limits(2));
        master_envelopes(freqband,elec_idx,:) = abs(hilbert(elec_bandpass));
    end
end

%% Save envelope mean trial matrices
for freqband = 1:length(FREQBAND)
            envelopes = squeeze(master_envelopes(freqband,:,:));
            env_filename = strcat(envelopeDir,'envelopes_',FREQBAND{freqband},'.mat');
            save(env_filename,'envelopes');
end

end
