function [gm_weight] = fn_gm_weight(primary_tissue,secondary_tissue)
%% Computes the portion of a bipolar electrode pair with GM signal
if strcmp(primary_tissue,'GM')
    if isempty(secondary_tissue)  % GM, GM
        gm_weight = 1;
    else                            % GM, WM; GM, OUT
        gm_weight = 0.6;
    end
else
    if isempty(secondary_tissue)  % WM, WM
        gm_weight = 0;
    else                            % WM, GM; OUT, GM
        gm_weight = 0.4;
    end
end

end