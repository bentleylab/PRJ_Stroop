
p_id = '2016-12-15_14-05-13';

py.path = ['C:\Users\KnightLab\Documents\Athina\StatReg\Python\', p_id, '\'];
p.path = ['C:\Users\KnightLab\Documents\Athina\StatReg\', p_id, '\'];
addpath('C:\Users\KnightLab\Documents\Athina\matlab\MatlabImportExport_v6.0.0\')
addpath('C:\Users\KnightLab\Documents\Athina\matlab\fieldtrip-20170207\')

s_id = 305;

cd(p.path)

xx = ls('trypy/*.ncs')
C = cellstr(xx);

ft_defaults
hdr = ft_read_header('trypy')
event = ft_read_event('trypy');

nSampleStart = 1;
nSampleStop = hdr.nSamples;
cfg                 = [];
cfg.dataset         = 'trypy';
cfg.trl             = [nSampleStart,nSampleStop,0];
data                = ft_preprocessing(cfg);

for kk = 1:length(data.hdr.chantype)
    
    data.hdr.chantype{kk} = 'elec';
    data.hdr.chanunit{kk} = '';
end

cfg=[];
cfg.resamplefs = 1000;
cfg.detrend    = 'no';
[data]        = ft_resampledata(cfg, data)
data.hdr.Fs = 1000;

%now we need to clean the triggers  - and only keep 10 or 5:
ev_clean = [];
for ev = 1:length(event)
    if event(ev).value == 10 |  event(ev).value == 5
        
        event(ev).sample = event(ev).sample / 8;
        ev_clean = [ev_clean event(ev)];
       
    end
end

data.cfg.event = ev_clean(1:end);

fiff_file  = [py.path, 'fif_P_',num2str(s_id), '.fif'];

fieldtrip2fiff(fiff_file, data)


