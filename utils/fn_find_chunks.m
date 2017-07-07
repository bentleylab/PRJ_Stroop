function [chunk_lim] = fn_find_chunks(data)
% Takes a binary (0/1) row vector and returns the edges of any chunks of data
% OUTPUTS:
%   chunk_lim [Nx2 array] - (start_ix,stop_ix) pairs for each chunk of consecutive data

chunk_edges = find(diff(data));
chunk_edges = [1 chunk_edges];
chunk_lim = NaN([numel(chunk_edges) 2]);
for chunk_ix = 1:length(chunk_edges)
    if chunk_ix==1
        chunk_lim(chunk_ix,1) = chunk_edges(chunk_ix);
    else
        chunk_lim(chunk_ix,1) = chunk_edges(chunk_ix)+1;  %+1 accoutns for lost sample when doing diff
    end
    if chunk_ix==length(chunk_edges)
        chunk_lim(chunk_ix,2) = size(data,2);
    else
        chunk_lim(chunk_ix,2) = chunk_edges(chunk_ix+1);
    end
end

end