function [i0refs, Vs, Y] = preprocess_current(data)
    [Vs, ~, ~] = unique(data{:,1});     % Get voltages
    [i0refs, ind, ~] = unique(data{:,2});
    nV = length(Vs);                         % Number of recorded voltages
    ni0 = length(i0refs);
    % Reshape so that rows = i_0,ref,f and columns = V
    y = data{:,end};
    Y = reshape(y, nV, ni0);
    Y = permute(Y, [2,1]);        % Swap columns to preserve index order
end