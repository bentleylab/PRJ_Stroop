function colors = fn_venn_colors(n_categories)
%% function colors = fn_venn_colors(n_categories)
% returns [n_cat, n_cat] matrix of colors for the number of categories
%   n=3 colors based on R, G, B (Y, P, C)
%   n=2 colors based on G, B (C)
% Stat color matrix:
%   Primary: Red, Green, Blue
%   Combinations: Yellow, Pink, Cyan
%   R   Y   P
%   Y   G   C
%   P   C   B

%% Define Colors
r = {[255 158 163]./255};
g = {[163 255 158]./255};
b = {[158 163 255]./255};
y = {[255 250 146]./255};
p = {[250 146 255]./255};
c = {[146 255 250]./255};

%% Select Colors
if n_categories==2
    colors = [g c; c b];
elseif n_categories==3
    colors = [r y p; y g c; p c b];
else
    error('n_categories must be 2 or 3');
end

end