function bp_lim = fn_bp_lim(fband)%SBJ,task,fband)
% Returns the upper and lower bounds for a given frequency band
% INPUTS: fband [str] = name of the frequency band of interest
% OUTPUTS: bp_lim [tuple int] = (lower_lim, higher_lim)

bp_lim = NaN([2 1]);
switch fband
    case 'delta'
        bp_lim(1) = 1; bp_lim(2) = 4;
    case 'thetaL'
        bp_lim(1) = 3; bp_lim(2) = 7;
        disp('WARNING: Are you sure you want thetaL (3-7 Hz) frequency band?');
    case 'theta'
        bp_lim(1) = 4; bp_lim(2) = 8;
    case 'theta-alpha'
        bp_lim(1) = 6; bp_lim(2) = 10;
    case 'alpha'
        bp_lim(1) = 8; bp_lim(2)= 12;
    case 'betaL'
        bp_lim(1) = 12; bp_lim(2)= 20;
    case 'betaH'
        bp_lim(1) = 20; bp_lim(2)= 35;
    case 'beta'
        bp_lim(1) = 15; bp_lim(2) = 30;
    case 'HG'
        bp_lim(1) = 70; bp_lim(2) = 150;
%     case {'cstULw','cstLw','cstHar','cstBeta'}
%         [bp_lim(1),bp_lim(2)]=func_returnCstFBand(SBJ,task,freqband);
    otherwise
        disp(strcat('Invalid frequency band variable: ',fband));
        error(strcat('Invalid frequency band variable: ',fband));
end

end
