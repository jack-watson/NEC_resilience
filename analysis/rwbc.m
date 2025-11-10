function [bc, info] = rwbc(G, W, useParallel)
% RWBC: Ridership-weighted betweenness centrality (nodes).

% [bc, info] = rwbc_passenger_betweenness(G, W)
% [bc, info] = rwbc_passenger_betweenness(G, W, useParallel)

% INPUTS
%   G           : MATLAB undirected graph. Uses:
%                 - G.Edges.Weight : nonnegative edge costs (distance).
%   W           : n-by-n OD passenger matrix (ordered pairs). Diagonal ignored.
%                 Sparse OK. Units arbitrary (passengers over the analysis period).
%   useParallel : (optional, logical) default false. If true and Parallel
%                 Toolbox is available, runs sources in parallel (parfor).

% OUTPUTS
%   bc : n-by-1 normalized passenger-weighted betweenness centrality.
%          Value is the fraction of all passengers whose shortest paths
%          traverse each node (excluding use as origin or destination).
%   info : struct with diagnostics:
%          .totalPassengers, .orderedPairs=true, .normalized=true,
%          .usedParallel, .tol, .notes

% METHOD
%   Deterministic geodesic RWBC:
%     C(v) = sum_{s!=v!=t} W(s,t) * (sigma_svt / sigma_st),
%   where sigma_st is the count of s->t shortest paths (by edge Weight),
%   and sigma_svt counts those paths that pass through v.
%   Implemented via a weighted Brandes accumulation with per-OD weights.

% REQUIREMENTS / ASSUMPTIONS
%   - G must be undirected with nonnegative weights.
%   - W size must match numnodes(G).
%   - Unreachable OD pairs contribute zero by construction.

% PERFORMANCE
%   - Prebuilds adjacency lists; per-source run is Dijkstra + Brandes accumulation.
%   - With ~1–2k nodes, parfor often helps if available.

% Core functionality/algorithm design hand-coded (JRW), ChatGPT-5 Pro used
% for optimization of time complexity (hence the style)

if nargin < 3 || isempty(useParallel), useParallel = false; end
assert(isa(G,'graph') && ~issymmetric(adjacency(G)) || isa(G,'graph'), ...
    'G must be a MATLAB graph.');
assert(ismember('Weight', G.Edges.Properties.VariableNames), ...
    'G.Edges.Weight is required (nonnegative costs).');

n = numnodes(G);
% Validate W
assert(all(size(W) == [n, n]), 'W must be n-by-n where n = numnodes(G).');
if ~issparse(W), W = sparse(double(W)); else, W = double(W); end
% Ignore self-trips
W = W - spdiags(diag(W), 0, n, n);

totalPassengers = full(sum(W(:)));
if totalPassengers <= 0
    bc = zeros(n,1);
    info = struct('totalPassengers',0,'orderedPairs',true,'normalized',true, ...
        'usedParallel',false,'tol',NaN,'notes',"Empty OD; returning zeros.");
    return
end

% Extract weighted adjacency (supports named nodes transparently).
ends = G.Edges.EndNodes;
if iscell(ends) || isstring(ends)
    % Map names to indices, if needed
    assert(ismember('Name', G.Nodes.Properties.VariableNames), ...
        'G.Nodes.Name is required when EndNodes are names.');
    [I,J] = mapEndNodesToIdx(G, ends);
else
    I = ends(:,1); J = ends(:,2);
end
w = double(G.Edges.Weight);
assert(all(isfinite(w) & w >= 0), 'Edge weights must be finite and >= 0.');

% Build symmetric adjacency lists (undirected graph).
nbrs = cell(n,1); wts = cell(n,1);
% Accumulate neighbors for both directions
deg = accumarray([I; J], 1, [n,1]);
deg = deg + accumarray([J; I], 1, [n,1]); %#ok<NASGU> % not directly used
% Preallocate via dynamic growth is OK for m-scale; keep it simple & robust.
for e = 1:numel(w)
    u = I(e); v = J(e); c = w(e);
    nbrs{u} = [nbrs{u}, v]; wts{u} = [wts{u}, c];
    nbrs{v} = [nbrs{v}, u]; wts{v} = [wts{v}, c];
end

% Distance equality tolerance for path counting (handle FP ties).
maxW = max(1, max(w));
tol  = 1e-12 * maxW * max(1, n);  % conservative scale with path length

% Accumulate RWBC
usedParallel = false;
if useParallel && license('test','Distrib_Computing_Toolbox')
    usedParallel = true;
    % Store per-source contributions, then reduce
    C_parts = zeros(n, n);   % n x n (OK for n ~ 2000; ~32 MB)
    parfor s = 1:n
        omega = full(W(s,:)).';  % target weights for this source
        C_parts(:,s) = brandes_source_weighted(nbrs, wts, s, omega, tol);
    end
    C = sum(C_parts, 2);
else
    C = zeros(n,1);
    for s = 1:n
        omega = full(W(s,:)).';
        C = C + brandes_source_weighted(nbrs, wts, s, omega, tol);
    end
end

% Normalize to fraction of all passengers
bc = C / totalPassengers;

info = struct();
info.totalPassengers = totalPassengers;
info.orderedPairs    = true;
info.normalized      = true;
info.usedParallel    = usedParallel;
info.tol             = tol;
info.notes           = "Deterministic geodesic; equal split across all shortest paths.";
end


% --------- helper: Brandes accumulation for one source with OD weights ---
function C_s = brandes_source_weighted(nbrs, wts, s, omega, tol)
    % omega(t) = W(s,t). Self-weight omega(s)=ignored automatically.
    n = numel(nbrs);
    % Dijkstra with predecessor tracking and sigma counts
    infv  = inf;
    dist  = infv(ones(n,1)); dist(s) = 0;
    sigma = zeros(n,1);      sigma(s) = 1;
    P     = cell(n,1);
    Sord  = zeros(n,1); top = 0; % stack of settled nodes (order of finalization)

    % Binary heap for Dijkstra
    heap.pos  = zeros(n,1);     % position in heap (0 = not in heap)
    heap.key  = zeros(n,1);
    heap.node = zeros(n,1);
    heap.size = 0;
    heap = heap_insert(heap, s, 0);

    while heap.size > 0
        [heap, u] = heap_pop(heap);            % settle u (min dist)
        top = top + 1; Sord(top) = u;

        nu = nbrs{u}; wu = wts{u};
        du = dist(u);

        for k = 1:numel(nu)
            v = nu(k);  alt = du + wu(k);
            dv = dist(v);
            if alt + tol < dv
                dist(v)  = alt;
                sigma(v) = sigma(u);
                P{v}     = u;                  % reset predecessors
                if heap.pos(v) == 0
                    heap = heap_insert(heap, v, alt);
                else
                    heap = heap_decrease_key(heap, v, alt);
                end
            elseif abs(alt - dv) <= tol
                sigma(v) = sigma(v) + sigma(u);
                P{v}     = [P{v}, u];
            end
        end
    end

    % Accumulation (Brandes): replace "1" by omega(w) for pair-weighting
    delta = zeros(n,1);
    C_s   = zeros(n,1);

    for ii = top:-1:1
        w = Sord(ii);
        % Skip unreachable nodes (dist = Inf) — they never get settled
        if sigma(w) == 0, continue; end

        coeff = (omega(w) + delta(w)) / sigma(w);
        % push dependencies to predecessors
        preds = P{w};
        for k = 1:numel(preds)
            v = preds(k);
            delta(v) = delta(v) + sigma(v) * coeff;
        end
        if w ~= s
            C_s(w) = C_s(w) + delta(w);  % endpoints excluded automatically: for t, delta(t)=0
        end
    end
end


% ---------------------------- binary heap utils --------------------------
function heap = heap_insert(heap, node, key)
    heap.size = heap.size + 1;
    i = heap.size;
    heap.node(i) = node;
    heap.key(i)  = key;
    heap.pos(node) = i;
    heap = sift_up(heap, i);
end

function heap = heap_decrease_key(heap, node, newkey)
    i = heap.pos(node);
    if newkey < heap.key(i)
        heap.key(i) = newkey;
        heap = sift_up(heap, i);
    end
end

function [heap, node] = heap_pop(heap)
    node = heap.node(1);
    last = heap.size;
    heap.pos(node) = 0;

    if last == 1
        heap.size = 0;
        return
    end
    heap.node(1) = heap.node(last);
    heap.key(1)  = heap.key(last);
    heap.pos(heap.node(1)) = 1;

    heap.size = last - 1;
    heap = sift_down(heap, 1);
end

function heap = sift_up(heap, i)
    while i > 1
        p = floor(i/2);
        if heap.key(i) < heap.key(p)
            [heap.node(i), heap.node(p)] = deal(heap.node(p), heap.node(i));
            [heap.key(i),  heap.key(p)]  = deal(heap.key(p),  heap.key(i));
            heap.pos(heap.node(i)) = i;
            heap.pos(heap.node(p)) = p;
            i = p;
        else
            break
        end
    end
end

function heap = sift_down(heap, i)
    while true
        l = 2*i; r = l + 1; m = i;
        if l <= heap.size && heap.key(l) < heap.key(m), m = l; end
        if r <= heap.size && heap.key(r) < heap.key(m), m = r; end
        if m ~= i
            [heap.node(i), heap.node(m)] = deal(heap.node(m), heap.node(i));
            [heap.key(i),  heap.key(m)]  = deal(heap.key(m),  heap.key(i));
            heap.pos(heap.node(i)) = i;
            heap.pos(heap.node(m)) = m;
            i = m;
        else
            break
        end
    end
end


% ----------- map EndNodes names to numeric indices (if needed) -----------
function [I,J] = mapEndNodesToIdx(G, ends)
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
end