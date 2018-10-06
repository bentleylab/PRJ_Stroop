
path_gen = '/Users/atzovara/Documents/Projects/StatReg/Irvine_data/';
patients = {'I875'};
file_ext = {'_0003'};

pp = 1;

p_id = char(patients(pp));

Timing.base = 100; % in ms
Timing.post = 500;
        
Timing.bin_l = 50; % bins in ms
Timing.base_bin = round(Timing.base/Timing.bin_l); % post-binned baseline length

micro.Names = {'mrth6'; 'mrth7'; 'mlth8'};
micro.clusters2test = [2,3,1]; % cluster we want to plot within each wire
micro.file_ext = char(file_ext{pp});
micro.p_id = p_id


path_spikes = [path_gen, p_id, filesep, 'Spikes', filesep];
events_file = [path_gen, p_id, filesep, 'events_Spikes.mat'];

Compute_FiringRates(path_spikes, events_file, micro,Timing)