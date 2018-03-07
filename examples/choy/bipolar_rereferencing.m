%% Rerefrence
%   data is a fiedltrip structure with the channels you want to reference
left_out_ch = {};   % This will report channels lost because they don't have a neighbor
for d = 1:numel(SBJ_vars.ch_lab.probes) %list of depth probes, e.g., {'RAM','ROF','LAC',...}
    cfg = [];
    cfg.channel = ft_channelselection(strcat(SBJ_vars.ch_lab.probes{d},'*'), data.label);
    probe_data = ft_selectdata(cfg,data);   % Grab data from this probe to plot in PSD comparison
    
    % Create referencing scheme
    if strcmp(SBJ_vars.ref_types{d},'BP')   %SBJ_vars.ref_types = cell with strings 'BP' or 'CAR' for each probe name
        cfg.montage.labelold = cfg.channel;
        [cfg.montage.labelnew, cfg.montage.tra, left_out_ch{d}] = fn_create_ref_scheme_bipolar(cfg.channel);
    elseif strcmp(SBJ_vars.ref_types{d},'CAR')
        cfg.reref      = 'yes';
        cfg.refchannel = setdiff(probe_data.label,SBJ_vars.ref_exclude);
        cfg.refmethod  = 'avg';
        left_out_ch{d} = {};    % CAR is applied to all channels
    else
        error(strcat('ERROR: Unrecognized reference type ',SBJ_vars.ref_types{d},...
            ' for probe ',SBJ_vars.ch_lab.probes{d}));
    end
    data_reref{d} = ft_preprocessing(cfg, data);
    
    if strcmp(psd_reref,'yes')  % Plot the PSD before and after for each channel, which nicely demonstrates the denoising
        psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/bp.reref/');
        if ~exist(psd_dir,'dir')
            mkdir(psd_dir);
        end
        fn_plot_PSD_1by1_compare_save(probe_data.trial{1},data_reref{d}.trial{1},...
            probe_data.label,data_reref{d}.label,data_reref{d}.fsample,...
            strcat(psd_dir,SBJ,'_PSD_bp.reref'),'bp','bp.reref',psd_fig_type);  %'bp'=bandpassed; 'bp.reref'=bandpass+reref
    end
end
clear data_bp probe_data

%Concatenate together again
cfg = [];
% cfg.appendsens = 'yes';
data = ft_appenddata(cfg,data_reref{:});
% Somehow, data.fsample is lost in certain cases here (new ft version I think):
data.fsample = data_reref{1}.fsample;
data_reref = data;

% Print left out channels, add to SBJ_vars
if ~isempty([left_out_ch{:}])
    fprintf('=============================================================\n');
    fprintf('WARNING: %i channels left out!\n',numel([left_out_ch{:}]));
    left_out_ch{:}
    fprintf('Consider adding these to SBJ_vars!\n');
    fprintf('=============================================================\n');
end
% SBJ_vars.ch_lab.left_out = [left_out_ch{:}];
