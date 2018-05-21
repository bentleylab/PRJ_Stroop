function data = rh_epoch(inputdata, trig, targettrig, response, blockno, fs, pre_t, post_t, label)

% function to epoch ECoG Data

% inputdata needs to be channel x samples
% trig is trigger x 2 with code and samplepoint
% targettrig is trigger to epoch on
% blockno is repition number

% 1. get onsets and trialno
onsets = trig(find(trig(:,1) == targettrig),2);
trialno = length(onsets);

% 2. create basic struct for FT
data.trial = {};
data.fsample = fs;
data.time = {};
data.label = label;
data.elec = [];
data.sampleinfo = zeros(trialno,2);
data.trialinfo = [ones(trialno,1)*targettrig, ones(trialno,1)*blockno, zeros(trialno,1)];
data.cfg = [];

% 3. loop over all trials and extract data
for i = 1:trialno
    
    startp  = onsets(i) - (pre_t * fs);
    endp    = onsets(i) + (post_t * fs);
    data.trial{1,i}(:,:) = inputdata(:,startp:endp);
    
    data.time{1,i} = linspace((-1*pre_t), post_t, (endp-startp+1));
    data.sampleinfo(i,[1 2]) = [startp endp];  
    
    % add reaction time from trig file
    targetpos = find(trig(:,2) == onsets(i));
    if trig(targetpos+1,1) == 7
        data.trialinfo(i,3) = (trig(targetpos+1,2) - trig(targetpos,2)) / fs; 
    else
        data.trialinfo(i,3) = NaN;
    end
    
end % end loop over trials

end