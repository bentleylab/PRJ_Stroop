% pop_resample() - resample dataset (pop up window).
%
% Usage:
%   >> [OUTEEG] = pop_resample( INEEG ); % pop up interactive window
%   >> [OUTEEG] = pop_resample( INEEG, freq);
%
% Graphical interface:
%   The edit box entitled "New sampling rate" contains the frequency in
%   Hz for resampling the data. Entering a value in this box  is the same 
%   as providing it in the 'freq' input from the command line.
%
% Inputs:
%   INEEG      - input dataset
%   freq       - frequency to resample (Hz)  
%
% Outputs:
%   OUTEEG     - output dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% Note: uses the resample() function from the signal processing toolbox
%       if present. Otherwise use griddata interpolation method (it should be
%       reprogrammed using spline interpolation for speed up).
%
% See also: resample(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: pop_resample.m,v $
% Revision 1.24  2007/02/20 14:18:37  arno
%  now can process several datasets
%
% Revision 1.23  2006/11/17 22:07:40  arno
% help message
%
% Revision 1.22  2006/11/17 21:58:07  arno
% fixed spline interpolation
%
% Revision 1.21  2006/11/17 21:51:18  arno
% *** empty log message ***
%
% Revision 1.20  2006/11/17 21:48:47  arno
% fixed resampling if sigproc absent
%
% Revision 1.19  2006/11/17 21:12:52  arno
% *** empty log message ***
%
% Revision 1.18  2005/10/01 23:09:51  arno
% remove doublet boundaries if any
%
% Revision 1.17  2005/09/27 22:09:30  arno
% now uses spline interpolation if the signal processing toolbox is absent (should be as efficient as the signal processing toolbox)
%
% Revision 1.16  2005/05/24 17:26:56  arno
% remove cell2mat
%
% Revision 1.15  2004/09/23 18:19:37  hilit
% corrected a bug, that in the case of boundaries in continuous data prevented from resampling.
%
% Revision 1.14  2004/09/15 06:03:24  arno
% now systematically crashes under Matlab 7?
%
% Revision 1.13  2004/08/10 22:54:52  arno
% fixed resampling for data epochs
%
% Revision 1.12  2004/08/03 01:33:03  arno
% convert to double for Matlab 7
%
% Revision 1.11  2004/06/03 19:18:18  arno
% taking into account boundaries for resampling
%
% Revision 1.10  2003/12/04 23:03:34  arno
% detect resample
%
% Revision 1.9  2003/12/03 03:28:34  arno
% message
%
% Revision 1.8  2003/12/03 03:19:23  arno
% allow to use interpolation
% instead of resample
%
% Revision 1.7  2003/09/29 23:45:07  arno
% reinitializing urevent if crash
%
% Revision 1.6  2003/06/28 02:28:05  arno
% fixing slight inacuracy in sampling rate
%
% Revision 1.5  2003/06/28 02:26:45  arno
% fixing slight inacuracy in new sampling rate
%
% Revision 1.4  2003/06/13 15:07:00  arno
% adding urevent compatiblity
%
% Revision 1.3  2003/02/16 22:59:40  arno
% adding text for GUI in header
%
% Revision 1.2  2002/08/12 02:30:30  arno
% [6~[6~inputdlg2
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-08-02 remove ica activity resampling (now set to []) -ad
% 03-08-02 debug call to function help -ad
% 04-05-02 recompute event latencies -ad

function [EEG, command] = pop_resample( EEG, freq); 

command = '';
if nargin < 1
    help pop_resample;
    return;
end;     
if isempty(EEG(1).data)
    disp('Pop_resample error: cannot resample empty dataset'); return;
end;    

if nargin < 2 

	% popup window parameters
	% -----------------------
	promptstr    = {['New sampling rate']};
	inistr       = { num2str(EEG(1).srate) };
	result       = inputdlg2( promptstr, 'Resample current dataset -- pop_resample()', 1,  inistr, 'pop_resample');
	if length(result) == 0 return; end;
	freq         = eval( result{1} );

end;

% process multiple datasets
% -------------------------
if length(EEG) > 1
    [ EEG command ] = eeg_eval( 'pop_resample', EEG, 'warning', 'on', 'params', { freq } );
    return;
end;

% finding the best ratio
[p,q] = rat(freq/EEG.srate, 0.0001); % not used right now 

% set variable
% ------------
EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
oldpnts   = EEG.pnts;

% resample for multiple channels
% -------------------------
if isfield(EEG, 'event') & isfield(EEG.event, 'type') & isstr(EEG.event(1).type)
    bounds = strmatch('boundary', { EEG.event.type });
    if ~isempty(bounds),
        disp('Data break detected and taken into account for resampling');
        bounds = [ EEG.event(bounds).latency ];
        if bounds(1) < 0, bounds(1) = []; end; % remove initial boundary if any
    end;
    bounds = [1 round(bounds-0.5)+1 size(EEG.data,2)+1];
    bounds(find(bounds(2:end)-bounds(1:end-1)==0))=[]; % remove doublet boundary if any
else 
    bounds = [1 size(EEG.data,2)+1]; % [1:size(EEG.data,2):size(EEG.data,2)*size(EEG.data,3)+1];
end;
if exist('resample') == 2
     usesigproc = 1;
else usesigproc = 0;
    disp('Signal Processing Toolbox absent: using custom interpolation instead of resample() function.');
    disp('This method uses cubic spline interpolation after anti-aliasing (see >> help spline)');    
end;

fprintf('resampling data %3.4f Hz\n', EEG.srate*p/q);
for index1 = 1:size(EEG.data,1)
    fprintf('.', index1);	
    sigtmp = reshape(EEG.data(index1,:, :), oldpnts, EEG.trials);
    
    if index1 == 1
        tmpres = [];
        indices = [1];
        for ind = 1:length(bounds)-1
            tmpres  = [ tmpres; myresample( double( sigtmp(bounds(ind):bounds(ind+1)-1,:)), p, q, usesigproc ) ];
            indices = [ indices size(tmpres,1)+1 ];
        end;
        if size(tmpres,1) == 1, EEG.pnts  = size(tmpres,2);
        else                    EEG.pnts  = size(tmpres,1);
        end;
        tmpeeglab = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
    else
        for ind = 1:length(bounds)-1
            tmpres(indices(ind):indices(ind+1)-1,:) = myresample( double( sigtmp(bounds(ind):bounds(ind+1)-1,:) ), p, q, usesigproc );
        end;
    end; 
    tmpeeglab(index1,:, :) = tmpres;
end;
fprintf('\n');	
EEG.srate   = EEG.srate*p/q;
EEG.data = tmpeeglab;

% recompute all event latencies
% -----------------------------
if isfield(EEG.event, 'latency')
    fprintf('resampling event latencies...\n');
    for index1 = 1:length(EEG.event)
        EEG.event(index1).latency = EEG.event(index1).latency * EEG.pnts /oldpnts;
    end;
    if isfield(EEG, 'urevent') & isfield(EEG.urevent, 'latency')
        try,
            for index1 = 1:length(EEG.event)
                EEG.urevent(index1).latency = EEG.urevent(index1).latency * EEG.pnts /oldpnts;
            end;
        catch, 
            disp('pop_resample warning: ''urevent'' problem, reinitializing urevents');
            EEG = rmfield(EEG, 'urevent');
        end;
    end;
end;

% resample for multiple channels ica
EEG.icaact = [];

% store dataset
fprintf('resampling finished\n');

EEG.setname = [EEG.setname ' resampled'];
EEG.pnts    = size(EEG.data,2);

command = sprintf('EEG = pop_resample( %s, %d);', inputname(1), freq);
return;

% resample if resample is not present
% -----------------------------------
function tmpeeglab = myresample(data, pnts, new_pnts, usesigproc)
    
    if usesigproc
        tmpeeglab = resample(data, pnts, new_pnts);
        return;
    end;
    
    % anti-alias filter
    % -----------------
    data         = eegfiltfft(data', 256, 0, 128*pnts/new_pnts);
    
    % spline interpolation
    % --------------------
    X            = [1:length(data)];
    nbnewpoints  = length(data)*pnts/new_pnts;
    nbnewpoints2 = ceil(nbnewpoints);
    lastpointval = length(data)/nbnewpoints*nbnewpoints2;        
    XX = linspace( 1, lastpointval, nbnewpoints2);
    
    cs = spline( X, data);
    tmpeeglab = ppval(cs, XX)';
