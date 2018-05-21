%% IR52 Processing Parameters
% Save out processing variables as proc_var after looking at cleaning output
SBJ = 'IR52';
import_dir = ['/home/knight/hoycw/PRJ_Stroop/data/' SBJ '/01_import/'];

%% Basics
preproc_vars.SBJ = SBJ;
preproc_vars.raw_filename = '2017010513_0022.besa';

%% Channel Selection
load(strcat(import_dir,SBJ,'_raw_labels.mat'));
preproc_vars.orig_n_ch = length(raw_labels);

preproc_vars.ch_lab.probes = {'RAM','RHH','RTH','RAC','ROF','RSMA','RAIN','RPIN','RBT',...
                         'LAM','LHH','LTH','LAC','LOF'};
preproc_vars.ch_lab.bad = {...
    'RAM3','RHH1','RHH2','RTH4','RTH5',...% epileptic
    'LAM3','LAM4','LHH3','LHH4','LTH3',...% epileptic
    'RAM1','RHH10','RTH1','RTH2','RAC10','LAM10','RAIN10',...% out of brain
    'RHH6',...%noisy
    'DC01','DC04'....% Not real data
    'REF','LUE','RUE','RLE',...% Not real data
    'EKG'...
    };
preproc_vars.ch_lab.eeg = {'FZ' 'CZ' 'OZ' 'C3' 'C4' 'O1' 'O2' 'T3' 'T4' 'F3' 'F4'};
preproc_vars.ch_lab.evnt = {...
    'DC02',...%photodiode
    'DC03',...%microphone
    };
preproc_vars.ch_lab.mic = {'DC03'};

%% Analysis Parameters
% most have only regular harmonics with normal width, and 120 and 240 are weak
% RBT, RPIN, RSMA have an extra peak at 200
% RHH6 has really bad at all harmonics, like LUE and other nonsense chan
preproc_vars.notch_freqs = [60 120 180 240 300]; %200 shoudl come out in re-referencing
preproc_vars.bs_width    = 2;

% Only 4 blocks, starts very quickly and ~75 extra seconds on the end
preproc_vars.analysis_time = {[0 510]};

%% Save proc_vars
out_filename = [import_dir SBJ '_preproc_vars.mat'];
save(out_filename,'preproc_vars');
