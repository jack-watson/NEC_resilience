function [T, edgeFlow, info] = gravityMdl(G, dCon, assignEdges, rshp)

%GRAVITYMDL  Gravity-model OD synthesis for an undirected rail graph

%   [T, edgeFlow, info] = gravityMdl(G)
%   [T, edgeFlow, info] = gravityMdl(G, dCon)
%   [T, edgeFlow, info] = gravityMdl(G, dCon, assignEdges)

% INPUTS
%   G            : MATLAB graph (undirected). Requires:
%                  - G.Nodes.rshp (vector of annual ons/offs per station)
%                  - G.Edges.Weight (edge distances/segment lengths)
%   dCon         : (optional, logical) if true, use doubly-constrained (Furness/IPF)
%                  Default: false (unconstrained; globally scaled)
%   assignEdges  : (optional, logical) if true, assign OD flows to edges via
%                  all-or-nothing (shortest paths) and return edgeFlow
%                  Default: false

% OUTPUTS
%   T        : n x n OD matrix (T(i,j): trips from i->j). diag(T)=0
%   edgeFlow : m x 1 edge flows (sum of OD flows along chosen shortest paths)
%              Empty if assignEdges=false. For undirected edges
%   info     : struct with metadata (alpha, model, balancing stats, etc.)

% MODEL
%   Power-law deterrence: f(d) = d^(-alpha) with alpha=1.5 (see below)
%   Productions = Attractions = G.Nodes.rshp
%   Unconstrained case scales a single K so sum(T)=sum(P)
%   Doubly-constrained uses classic IPF:
%     T_ij = a_i * b_j * P_i * A_j * f(d_ij), enforcing row/col sums

% NOTES
% NOTES
% - No self-trips (diag(T)=0)
% - Unreachable pairs (d=Inf) assigned 0 flow
% - No symmetry enforced (T_ij ~= T_ji allowed)
% - Edge assignment uses weighted shortest-path trees per origin for efficiency
% - OD build is vectorized; edge assignment is O(n * (n + m)) with small 
%      constants (tree-based aggregation per origin)

% Combination of hand-coding (JRW) and ChatGPT-5 Pro
% Jack R. Watson, 2025

% ------------------------ argument checking ------------------------
if nargin < 2 || isempty(dCon)
    dCon = false;
end
if nargin < 3 || isempty(assignEdges)
    assignEdges = false;
end

assert(isa(G,'graph'), 'G must be a MATLAB graph (undirected).');
if nargin < 4
    assert(ismember('rshp',   G.Nodes.Properties.VariableNames), ...
        'G.Nodes.rshp (ridership) is required.');
end
assert(ismember('Weight', G.Edges.Properties.VariableNames), ...
    'G.Edges.Weight (edge distance) is required.');

n = numnodes(G);
m = numedges(G);
if nargin < 4 || isempty(rshp)
    P = double(G.Nodes.rshp(:)); % productions
else
    P = double(rshp);
end
A = P; % attractions (as specified)
totalP = sum(P);
if totalP <= 0 || n == 0
    T = zeros(n); edgeFlow = zeros(m,1);
    info = struct('alpha',NaN,'model','gravity-power','dCon',dCon, ...
        'totalTrips',0,'K',NaN,'converged',true);
    return
end

% Deterrence exponent (power law). Adjust here if desired.
alpha = 1.5;

% ------------------------ distances & deterrence -------------------
% All-pairs shortest-path distances using edge weights.
% (distances(G) uses G.Edges.Weight if present)
D = distances(G);
% Build deterrence matrix F = d^(-alpha); zero for diag and Inf
F = zeros(n, n);
finiteMask = isfinite(D) & ~eye(n);
F(finiteMask) = D(finiteMask) .^ (-alpha);
% unreachable pairs => F = 0; diagonal => F = 0

% ------------------------ OD synthesis -----------------------------
if ~dCon
    % Unconstrained gravity with global scaling K so that sum(T)=sum(P)
    base = (P * A.');     % outer product
    base = base .* F;     % apply deterrence & zero diag/unreachables
    S = sum(base(:));     % denominator
    if S <= 0
        T = zeros(n);
        infoK = 0;
    else
        infoK = totalP / S;
        T = infoK * base;
    end
    converged = true; iters = 0;
else
    % Doubly constrained: Furness/IPF for a_i, b_j
    % T_ij = a_i * b_j * P_i * A_j * F_ij, s.t. rows sum to P, cols to A.
    a = ones(n,1);
    b = ones(n,1);
    tol = 1e-8;
    maxIter = 30000;
    converged = false;
    % Pre-handle zero P or A to avoid divisions
    posP = P > 0;
    posA = A > 0;

    % Iterative proportional fitting
    for iters = 1:maxIter
        % Update a:  a_i = 1 / sum_j (b_j * A_j * F_ij)
        denom_row = F * (b .* A);
        a_new = zeros(n,1);
        idx = (denom_row > 0) & posP;
        a_new(idx) = 1 ./ denom_row(idx);
        a = a_new;

        % Update b:  b_j = 1 / sum_i (a_i * P_i * F_ij)
        denom_col = F.' * (a .* P);
        b_new = zeros(n,1);
        idx = (denom_col > 0) & posA;
        b_new(idx) = 1 ./ denom_col(idx);
        b = b_new;

        % Construct T and check convergence
        T = (a .* P) * ( (b .* A).' );
        T = T .* F;  % enforce deterrence and zeros

        rowErr = max(abs(sum(T,2) - P)) / max(1, totalP);
        colErr = max(abs(sum(T,1).' - A)) / max(1, totalP);
        if max(rowErr, colErr) < tol
            converged = true;
            break
        end
    end
    infoK = NaN;
end

% Ensure strict zero diagonal
T(1:n+1:end) = 0;

% Diagnostics
info = struct();
info.alpha        = alpha;
info.model        = 'gravity-power';
info.dCon         = dCon;
info.K            = infoK;
info.totalTrips   = sum(T(:));
info.iters        = exist('iters','var') * iters;
info.converged    = converged;
info.rowL1Err     = norm(sum(T,2) - P, 1);
info.colL1Err     = norm(sum(T,1).' - A, 1);

if ~assignEdges % edge assignment
    edgeFlow = [];
    return
end

% Map undirected edge (u,v) -> edge index
ends = G.Edges.EndNodes; % m-by-2
if iscell(ends) || isstring(ends)
    % Require a Name column to map names -> indices
    assert(ismember('Name', G.Nodes.Properties.VariableNames), ...
        ['G.Nodes.Name is required to map EndNodes names to indices. ' ...
        'Add a unique Name for each node.']);

    nodeNames = string(G.Nodes.Name);
    endNames  = string(ends);

    [tfI, I] = ismember(endNames(:,1), nodeNames);
    [tfJ, J] = ismember(endNames(:,2), nodeNames);

    if ~all(tfI & tfJ)
        missing = unique([endNames(~tfI,1); endNames(~tfJ,2)]);
        error('EndNodes contain unknown names not found in G.Nodes.Name: %s', ...
            strjoin(cellstr(missing), ', '));
    end

    I = double(I); J = double(J);
else
    % Numeric indices
    I = ends(:,1); J = ends(:,2);
end

EIdx = sparse([I; J], [J; I], [(1:m)'; (1:m)'], n, n);

% Precompute distance tolerance for parent selection (to handle fp error)
Dfin = D(isfinite(D) & D>0);
if isempty(Dfin)
    distTol = 1e-12;
else
    distTol = 1e-12 * max(Dfin);
end

edgeFlow = zeros(m,1);

% For efficiency, reuse D for each origin s (single-source tree logic)
% Aggregate flows on shortest-path tree from s using subtree summation.
for s = 1:n
    d = D(s,:).';
    % Skip completely isolated origins
    if ~isfinite(d(s))
        continue
    end

    % Build parent for each node using distance optimality:
    % parent(v) is neighbor u of v s.t. d(u) + w(u,v) == d(v) (within tol),
    % with preference for neighbor closer to s.
    parent = zeros(n,1); parent(s) = s;
    for v = 1:n
        if v == s || ~isfinite(d(v)), continue; end
        nbrs = neighbors(G, v);
        if isempty(nbrs), continue; end
        eids = findedge(G, v*ones(numel(nbrs),1), nbrs);
        w    = G.Edges.Weight(eids);
        cand = abs(d(nbrs) + w - d(v)) <= distTol & (d(nbrs) < d(v) + distTol);
        if any(cand)
            % choose the neighbor with smallest d (closer to s)
            nbrs_c = nbrs(cand);
            [~,k]  = min(d(nbrs_c));
            parent(v) = nbrs_c(k);
        else
            % if no candidate (rare due to numeric issues), pick the best approx
            [~,k] = min(abs(d(nbrs) + w - d(v)));
            parent(v) = nbrs(k);
        end
    end

    % Subtree aggregation: start with flows from s to each destination v
    subflow = T(s,:).';
    subflow(~isfinite(d)) = 0; % unreachable
    subflow(s) = 0;

    % Process nodes in decreasing distance from s (leaves to root)
    [~, ord] = sort(d, 'descend');
    for kk = 1:n
        v = ord(kk);
        if v == s || parent(v) == 0 || ~isfinite(d(v)), continue; end
        f = subflow(v);
        if f == 0, continue; end
        eid = EIdx(v, parent(v));

        if eid ~= 0
            edgeFlow(eid) = edgeFlow(eid) + f;
        end

        subflow(parent(v)) = subflow(parent(v)) + f;
    end

end

end




