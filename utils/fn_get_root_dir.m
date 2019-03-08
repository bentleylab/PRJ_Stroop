function [root_dir, app_dir] = fn_get_root_dir()
%% Check with user/OS/root directory, return relevant locations
% OUTPUTS:
%   root_dir [str] - source dir with PRJ_??? in it
%   app_dir [str] - source of fieldtrip, other apps

if exist('/home/knight/hoycw/','dir')
    root_dir = '/home/knight/hoycw/';
    app_dir   = [root_dir 'Apps/'];
elseif exist('/Volumes/hoycw_clust/','dir')
    root_dir = '/Volumes/hoycw_clust/';
    app_dir   = '/Users/colinhoy/Code/Apps/';
else
    error('root directory not found. where are you running this?');
end

end