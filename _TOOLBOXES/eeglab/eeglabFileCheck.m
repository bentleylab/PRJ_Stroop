function [eegStruct] = eeglabFileCheck(fileExtension, eegDir)

%
% function [eegStruct] = eeglabFileCheck(fileExtension)
%
% This function requires one input: the three-character file extention of
% the files you want to look at within a directory (as a char array).
%
% This function provides one output, a simple array that assigns subject
% numbers to each file in the directory ('01', '02, '03'...) and the
% associated filename.

% Bradley Voytek
% Copyright (c) 2007
% University of California, Berkeley
% Helen Wills Neuroscience Institute
% btvoytek@berkeley.edu

%eegDir = [pwd '/'];

% Load file names into array
eegStruct = [];
eegFiles = dir(eegDir);
if(isempty(eegFiles))
  fprintf('\n\nERROR: No files found in eeglabFileCheck: [%s] [%s]\n\n', eegDir, fileExtension);
  return;
end

% Determines if files in a directory are *.<fileExtension> files and, if 
% so, adds a "process = true" tag to them
for eegIt = 1:length(eegFiles)
        if  eegFiles(eegIt).isdir == 0
            eegStr = eegFiles(eegIt).name((length(eegFiles(eegIt).name)-length(fileExtension)):(length(eegFiles(eegIt).name)));
            if strcmp(eegStr, ['.' fileExtension]) == 1
                eegFiles(eegIt).process = 1;
            elseif strcmp(eegStr, ['.' fileExtension]) == 0
                eegFiles(eegIt).process = 0;
            end
        elseif eegFiles(eegIt).isdir == 1
            eegFiles(eegIt).process = 0;
        end
end
clear eegIt eegStr;
% **********

% Create filename structure containing only *.<fileExtension> files
procCount = 1;

for eegIt = 1:length(eegFiles)
    if eegFiles(eegIt).process == 1
        eegStruct(procCount).name = eegFiles(eegIt).name;
        procCount = procCount + 1;
    end
end
clear eegIt eegFiles procCount;
% **********

% Add leading '0' to subject number for file-naming purposes
% for zeroIt = 1:length(eegStruct)
%     if zeroIt < 10
%         eegStruct(zeroIt).number = ['0' num2str(zeroIt)];
%     elseif zeroIt >= 10
%         eegStruct(zeroIt).number = [num2str(zeroIt)];
%     end
% end
% 
% clear zeroIt;
% **********