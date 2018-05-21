%% Filedtrip to MNE example code

% Export FT to MNE Raw
cfg = [];
cfg.dataset = 'MarkusBraille.ds';
cfg.trialdef.triallength = Inf;
cfg = ft_definetrial(cfg);
 
cfg.continuous = 'yes';
cfg.channel = {'MEG', '-MLP31', '-MLO12'};
data = ft_preprocessing(cfg);

fiff_file  = 'ctf_raw.fif';
fieldtrip2fiff(fiff_file, data)

% Import MNE raw to FT
fiff_file = 'mne_python_raw.fif';
 
cfg = []
cfg.dataset = fiff_file;
data1 = ft_preprocessing(cfg);
ft_datatype(data1)  % returns 'raw'
 
event = mne_read_events('ctf_raw-eve.fif')

