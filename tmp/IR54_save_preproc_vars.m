%% IR52 Pre-processing Parameters
% Save out processing variables as proc_var after looking at cleaning output
SBJ = 'IR54';
import_dir = ['/home/knight/hoycw/PRJ_Stroop/data/' SBJ '/01_import/'];

%% Basics
preproc_vars.SBJ = SBJ;
preproc_vars.raw_filename = '2017012411_0008.besa';

%% Channel Selection
load(strcat(import_dir,SBJ,'_raw_labels.mat'));
preproc_vars.orig_n_ch = length(raw_labels);

preproc_vars.ch_lab.probes = {'RAM','RHH','RTH','RAC','ROF',...
                                'LAM','LHH','LTH','LAC','LOF'};
preproc_vars.ch_lab.bad = {...
    'RAM1','RAM2','RAM8','RAM9','RAM10','RHH1','RHH2','RHH3','RHH4',...% epileptic
    'ROF1','LAM10',...% bad- both are loose
    'ROF10','LTH10','LAC10','LOF1',...% out of brain
    'RSMA*','RBT*','RAIN*','RPIN*',...% Not real data
    'DC01','DC04'....% Not real data
    'REF','LLE','LUE','RLE','RUE',...% Not real data
    'EKG'...
    };
%     'LAC5','LAC6','LAC7','LAC8','LOF5','LOF6','LOF7','LOF8','LOF9','LOF10',...
%     slowing noted in other datasets, but Bob says it's okay in this one
%     and I can keep them; that said, WATCH OUT!!!
preproc_vars.ch_lab.eeg = {'FZ' 'CZ' 'OZ' 'C3' 'C4'};
preproc_vars.ch_lab.evnt = {...
    'DC02',...%photodiode
    'DC03',...%microphone
    };
preproc_vars.ch_lab.mic = {'DC03'};

%% Analysis Parameters
% pretty standard
% 100 and 200 are small but consistent
preproc_vars.notch_freqs = [60 120 180 240 300]; %100 200 come out with rereferencing
preproc_vars.bs_width    = 2;

% Data starts ~107s, ends ~1101s; end of recording is ~1114s 
preproc_vars.analysis_time = {[90 1114]};

%% Save proc_vars
out_filename = [import_dir SBJ '_preproc_vars.mat'];
save(out_filename,'preproc_vars');




% % Remove extra characters from channel labels
% for channel_n = 1:length(data.label)
%     data.label{channel_n} = strrep(data.label{channel_n},'POL','');
%     data.label{channel_n} = strrep(data.label{channel_n},'Ref','');
% end

