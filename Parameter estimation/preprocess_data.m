function [X1,X2,Y] = preprocess_data(data)
    [Vs,    indV, ~] = unique(data{:,1});       % Get voltages
    [I_0a,  ind0a,~] = unique(data{:,2});       % Get all i_0a
    [I_0c,  ind0c,~] = unique(data{:,3});       % Get all i_0c
    nV = length(Vs);                         % Number of recorded voltages
    na = length(I_0a);                       % Number of i_0a
    nc = length(I_0c);                       % Number of i_0c
    % Create grid
    % 1st index = i_0a
    % 2nd index = i_0c
    [X1,X2] = ndgrid(I_0a,I_0c);
    % Reshape average current values to accomodate grid
    y = data{:,end};
    Y = reshape(y, nV, nc, na);
    Y = permute(Y, [1,3,2]);        % Swap columns to preserve index order
end