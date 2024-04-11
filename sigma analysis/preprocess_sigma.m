function [sigmas, Vs, Y] = preprocess_sigma(data)
    [Vs, ~, ~] = unique(data{:,1});     % Get voltages
    [sigmas, ind, ~] = unique(data{:,2});
    nV = length(Vs);                         % Number of recorded voltages
    nsigma = length(sigmas);
    % Reshape so that rows = s_io and columns = V
    y = data{:,end};
    Y = reshape(y, nV, nsigma);
    Y = permute(Y, [2,1]);        % Swap columns to preserve index order
end