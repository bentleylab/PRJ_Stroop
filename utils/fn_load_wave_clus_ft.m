function [spike] = fn_load_wave_clus_ft(SBJ,cluster_dir)
%% This function loads times_chan.mat outputs from wav_clus spiking and 
%   converts them into fieldtrip compatible spike structs
%   NOTE: This function takes the place of ft_read_spike because I want to
%   use the continuous data post-clustering, rather than the .nse files.
% INPUTS:
%   SBJ [str] - name of the patient (to load SBJ_vars, e.g., analysis_time)
%   cluster_dir [str] - full path to the clustering output files (e.g., times_channel.mat)
% OUPUTS:
%   spike [ft struct] - struct with the following fields:
%       spike.label     = 1xNchans cell-array, with channel labels
%       spike.waveform  = 1xNchans cell-array, each element contains a matrix (Nleads x Nsamples X Nspikes)
%       spike.waveformdimord = '{chan}_lead_time_spike'
%       spike.timestamp = 1xNchans cell-array, each element contains a vector (1 X Nspikes)

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath(ft_dir);
ft_defaults

%% Directories
% SBJ_vars
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
if numel(SBJ_vars.analysis_time)>1
    error('havent set up multi-run processing yet!');
elseif numel(SBJ_vars.analysis_time{1})>1
    error('havent set up processing for multi block concat!');
end

% Get file names
if strcmp(cluster_dir,'')
    cluster_dir = [SBJ_vars.dirs.preproc 'micro_clusters/semi_auto/'];
end
clust_fnames = dir([cluster_dir 'times*' SBJ_vars.ch_lab.suffix '.mat']);

% Initialize spike struct
spike.label = {};
spike.timestamp = {};
spike.waveform = {};
spike.waveformdimord = '{chan}_lead_time_spike';

%% Load cluster data
unit = 0;
for ch = 1:numel(clust_fnames)
    % Get Channel label
    uscores = strfind(clust_fnames(ch).name,'_');
    ch_lab  = clust_fnames(ch).name(uscores(1)+1:uscores(2)-1);
    % Get first time stamp
    if ch==1
        micro_hdr = ft_read_header([SBJ_vars.dirs.SU 'micro/' ch_lab SBJ_vars.ch_lab.suffix '.ncs']);
    end
    % Load clusters
    load([cluster_dir clust_fnames(ch).name]);
    
    %% Process Channel Clusters
    clust_n = sort(unique(cluster_class(:,1))>0);
    % Remove rejected cluster
    clust_n(clust_n==0) = [];
    
    % Adjust spike times from 1000 kHz time stamps to sec
    cluster_class(:,2) = (cluster_class(:,2)-double(micro_hdr.FirstTimeStamp)/1000)/1000;
    % Account for analysis time
    cluster_class(:,2) = cluster_class(:,2) - SBJ_vars.analysis_time{1}{1}(1);
    % Toss spikes before analysis_time
    cluster_class = cluster_class(cluster_class(:,2)>0,:);
    
    % Gather spike information for each cluster/unit
    for cl = 1:numel(clust_n)
        unit = unit+1;
        unit_idx = cluster_class(:,1)==cl;
        
        spike.label{unit} = [ch_lab '-' num2str(cl)];
        spike.timestamp{unit} = cluster_class(unit_idx,2)';
        spike.waveform{unit}  = zeros([1 size(spikes,2) sum(unit_idx) ]);    % Create singleton dimension to fit dimord
        spike.waveform{unit}(1,:,:)  = spikes(unit_idx,:)';
    end
end

end