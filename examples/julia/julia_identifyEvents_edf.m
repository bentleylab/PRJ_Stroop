function [trl, event] = identifyEvents_edf(cfg)
% function [trl, event] = identifyEvents_edf(cfg,triggerIndex,triggerDirection,threshold)

% Read the header information
hdr                 = ft_read_header(cfg.dataset);

% Set i) photodiode chan number & ii) threshold for signal changes in photodiode
chanindx        = 2; % channel number of photodiode
detectflank     = 'both'; % does signal go up to indicate event (white square on black bckgrd)
threshold       = 100000; %'3*nanmedian'; 100000 for IR
event           = ft_read_event(cfg.dataset, 'chanindx', chanindx, 'detectflank', detectflank, 'threshold', threshold);

% To figure out the best threshold for detecting photodiode signals...
% dbstop if error
% set debugging stopper at line 68 in read_trigger
% run ft_definetrials
% plot dat

trl = [];
waitfordown = 0;

pretone = 1 * hdr.Fs; % in secs
posttone = 2 * hdr.Fs; % in secs
preswitch = 1 * hdr.Fs; % in secs
postswitch = 4 * hdr.Fs; % in secs

for i = 1:length(event)
    
    if ~isempty(strfind(event(i).type, 'up')) % trigger onset

        goneup = i;
        waitfordown = 1;
        
    elseif ~isempty(strfind(event(i).type, 'down')) && waitfordown % trigger offset
        
        onset = event(goneup).sample;
        offset = event(i).sample;
        
        if (event(i).sample-event(goneup).sample) < (0.5 * hdr.Fs) % tone
            trlbegin   = event(goneup).sample - pretone;
            trlend     = event(goneup).sample + posttone - 1;
        else
            trlbegin   = event(goneup).sample - preswitch;
            trlend     = event(goneup).sample + postswitch - 1;
        end
        
        newtrl    = [trlbegin trlend onset offset];
        trl       = [trl; newtrl]; % store in the trl matrix
        
    end
    
end

