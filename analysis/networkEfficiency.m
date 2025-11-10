function E = networkEfficiency(y, params)
% NETWORKEFFICIENCY computes the global network efficiency at a given time,
% based on the underlying L-space (physical) graph.

% IN:
%   y: A vector of node integrities at a given time.
%   params: Structure with the following fields:
%           - G: A MATLAB graph object representing the physical infrastructure.
%           - numNodes: Total number of nodes in the network.

% The induced subgraph of functional nodes is obtained and the standard definition of global efficiency is applied:
%    E(G) = (1/(N*(N-1))) * sum_{i ≠ j} 1/d_{ij},
% where the sum is over all pairs of functional nodes (pairs including non-functional nodes contribute 0).

if nargin == 1
    params = struct;
    params.G = y;
    params.numNodes = height(params.G.Nodes);
    y = ones(height(params.G.Nodes),1);
    if all(params.G.Edges.Weight == mean(params.G.Edges.Weight))
        params.isWeighted = false;
    else
        params.isWeighted = true;
    end
    params.resOpts.weightByRshpTF = false;
end

N = params.numNodes; % total number of nodes in the full network
G = params.G;
if isfield(params, "isWeighted") && params.isWeighted
    isWeighted = true;
else
    isWeighted = false;
end

ydims = size(y);
if ydims(1) > 1 && ydims(2) > 1
    mask   = y==1; % logical once
    nCol   = size(mask,2);
    E     = zeros(1,nCol);
    for k = 1:nCol
        idx = mask(:,k);
        if k > 1
            if isequal(idx, lastIdx)
                E(k) = E(k-1);
                continue
            end
        end
        lastIdx = idx;
        Nk  = nnz(idx);
        if Nk < 2 % 0 or 1 nodes --> efficiency = 0
            E(k) = 0;
            continue
        end
        H = subgraph(G, idx); % induced sub‑graph
        if isWeighted
            D = distances(H); % auto‑selects Dijkstra for weights
        else
            D = distances(H,'method','unweighted'); % breadth-first search
        end
        % Convert to reciprocal, ignore self‑loops and Infs
        D(1: Nk+1 : end)  = inf; % discard zeros on diagonal
        invD              = 1./D; % Inf --> 0 automatically
        invD(invD == Inf) = 0;
        E(k)              = sum(invD(:)) / (N*(N-1)); % should always be considering N as the total number of nodes in G0, NOT the number of undamaged nodes at timestep k
    end
else

    F = find(y == 1); % Identify functional nodes
    if numel(F) < 2  % If fewer than two nodes are functional, set efficiency = 0
        E = 0;
        return;
    end

    subG = subgraph(params.G, F); % Induce subgraph of functional nodes
    D = distances(subG); % Compute all-pairs shortest-path distances on the induced graph
    D(isinf(D)) = 0; % Handle any infinite distances (disconnected pairs) by assigning 0


    invD = 1 ./ D; % Compute inverse distances
    invD(isinf(invD)) = 0;
    invD(eye(size(invD))==1) = 0; % Remove self-loops by ensuring the diagonal is zero

    if params.resOpts.weightByRshpTF
        % Always uses indv; G indv is set to G NoN in invoking function seqFailRecovNEC (for corridor only, not for indv systems ofc)
        ODraw = params.resOpts.ODMatrix.indv;
        ODsum = params.resOpts.ODMatrix.sum.indv;

        isActiveTF = logical(params.resOpts.ODMatrix.activeNodeIndices);
        ODraw(~isActiveTF,:) = [];
        ODraw(:,~isActiveTF) = [];
        ODraw(ODraw <= 0) = 1; % enforces no zero or negative ridership other than diagonal
        ODraw = ODraw.*double(eye(size(ODraw,1)) == 0); % enforce zeros on diagonal (no self-ridership)
        OD = rescale(ODraw ./ ODsum, 1, 2);
        Eij = invD.*OD;
    else
        Eij = invD;
    end

    sum_inv = sum(Eij(:)); % Sum contributions from all pairs
    E = sum_inv / (N*(N-1)); % Normalize by total number of ordered pairs in full network

end

end




