function G = get_communicability(adjacency,type)
% A: n*n adjacency matrix
% type: 0, binary (default); 1, weighted.
% G: n*n communicability matrix

if ~exist('type','var')
    if isempty(setdiff([0,1], unique(adjacency)))
        type = 0;
    else
        type = 1;
    end
end

switch type
    case 0
        G = expm(adjacency);
    case 1
        % negative square root of nodal degrees
        row_sum = sum(adjacency, 2);    % Calculate row sums of adjacency matrix
        neg_sqrt = row_sum .^ -0.5;     % Calculate the negative square root of nodal degrees
        square_sqrt = diag(neg_sqrt);   % Create a diagonal matrix with negative square roots as diagonal elements

        % normalize input matrix
        for_expm = square_sqrt * adjacency * square_sqrt;  % Normalize adjacency matrix

        % calculate matrix exponential of normalized matrix
        G = expm(for_expm);
end

G(logical(eye(size(G)))) = 0; % Set diagonal elements to zero
