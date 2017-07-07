function PLV_time = fn_PAC_PLVenvPh_ePair_time(phase1, phase2, starts, ends, varargin)%, foldername, elec1, elec2)
%Taking 2 electrodes and calculating average phase coherence over all trials
% This version can take the arg 'return_complex', which if 1 will return a
% complex time series instead of the absolute value of that time series

% OPTION: pass 2 extra args to get the complex time series (without taking abs()) back
%   e.g., func_PAC_PLVenvPh_ePair_time(phase1, phase2, starts, ends, 'return_complex', 1)
for n=1:2:length(varargin)-1
    switch lower(varargin{n})
        case 'return_complex'
            return_complex = varargin{n+1};
%             disp(strcat('func_PAC_PLVenvPh_ePair_time VARARGIN RECIEVED: return_complex = ',num2str(return_complex)));
    end
end
if ~exist('return_complex', 'var')
    return_complex = 0;
end

% Cutoff signals by time range
ch1 = NaN(length(starts),length(phase1(starts(1):ends(1)))); %all wins should be same length, so this should work
ch2 = NaN(length(starts),length(phase1(starts(1):ends(1))));
for k = 1:length(starts)
    a = starts(k);
    b = ends(k);
    ch1(k,:) = phase1(a:b);
    ch2(k,:) = phase2(a:b);
end

% Angle diff at each time point for each trial
angle_diff = ch1 - ch2;

%PLV_trials averages across time (get one PLV for each trial)
% PLV_vector = abs(mean(exp(1i*angle_diff),2));
% PLV_across_trials = mean(PLV_vector);

%PLV_time averages across trials (get a PLV time series with one PLV for each time point)
if return_complex
    PLV_time = mean(exp(1i*angle_diff),1);
else
    PLV_time = abs(mean(exp(1i*angle_diff),1));
end
% %Save PLV across trials
% cd(foldername);
% save(strcat('PAC_PLVenvPh_trials_e',num2str(elec1,'%03d'),'_e',num2str(elec2,'%03d')), 'PLV_vector');
% save(strcat('angleDiff_trialXtime_e',num2str(elec1,'%03d'),'_e',num2str(elec2,'%03d')), 'angle_diff');

end




