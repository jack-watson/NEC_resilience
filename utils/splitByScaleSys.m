function [S, paramsS, G] = splitByScaleSys(G, params)
%SPLITBYSCALESYSTEM  Split an undirected graph by single-value (scale,system)
% pairs, letting multi-membership nodes join the relevant sub-graphs.

%   Each entry in G.Nodes.scale / G.Nodes.system is either
%       – a single char/string     --> treated as that lone value
%       – a cell array of char/str --> treated as a list of values
%   A sub-graph is created *only* for (scale,system) pairs that appear
%     in at least one node whose attribute lists are BOTH length-1.
%   A node belongs to sub-graph (s,sys) if its scale list includes *s*
%     AND its system list includes *sys* (case-insensitive).
%   All edges among included nodes are retained.
%   Output S is a cell array of MATLAB graph objects; order arbitrary.

%% 1. Normalise attributes → cell array of lower-case string vectors
n            = numnodes(G);
nodeScales   = cell(n,1);
nodeSystems  = cell(n,1);

for i = 1:n
    nodeScales{i}  = toStringVec(G.Nodes.scale{i});
    nodeSystems{i} = toStringVec(G.Nodes.system{i});
end

%% 2. Collect “atomic” (scale,system) pairs from nodes with length-1 lists

atomicPairs = cell(n,2);
hasScalarPair = false(n,1);

for i = 1:n
    atomicPairs{i,1} = nodeScales{i};
    atomicPairs{i,2} = nodeSystems{i};
    if isscalar(nodeScales{i}) && isscalar(nodeSystems{i})
        hasScalarPair(i) = true;
    end
end

[pairs, ~, pidx] = unique(string(atomicPairs(hasScalarPair,:)), "rows");   % P-by-2 cell array of strings
if isempty(pairs)
    S = {};   % nothing to split
    return
elseif nargout >= 3
    G.Nodes.scaleSysPair = atomicPairs;
    G.Nodes.pairIndex = zeros(n,1);
    G.Nodes.pairIndex(hasScalarPair) = pidx;
end

%% List of bad edges to remove before returning processed subsystems
% 1. CT Rail
e2rmCTR = ["Bridgeport", "New Haven"; "Old Saybrook", "New Haven"; "Bridgeport", "Stamford"];

% 2. LIRR - ALL GOOD!

% 3. MARC
e2rmMARC = ["Martinsburg", "Harpers Ferry"; "Harpers Ferry", "Rockville"; "Union Station, Washington, DC", "Rockville";...
    "Penn Station, Baltimore", "BWI Airport"; "BWI Airport", "New Carrollton"; "Penn Station, Baltimore", "Aberdeen"];

% 4. MBTA Commuter - ALL GOOD!

% 5. Metro-North
e2rmMNR = ["Croton-Harmon", "Poughkeepsie"; "Croton-Harmon", "Yonkers"; "New Rochelle", "Stamford"; "Stamford", "South Norwalk";...
    "Stamford", "Bridgeport"; "South Norwalk", "Bridgeport"; "Bridgeport", "New Haven"];

% 6. NJ Transit
e2rmNJT = ["Trenton", "Princeton Junction"; "Princeton Junction", "New Brunswick"; "New Brunswick", "Metropark";...
    "Metropark", "Liberty International"; "Penn Station Newark", "Penn Station NYC"];

% 7. SEPTA Regional
e2rmSEPR = ["Trenton", "Cornwells Heights"; "Cornwells Heights", "North Philadelphia"; "Trenton", "Princeton Junction";...
    "30th St", "Wilmington"; "Wilmington", "Newark, DE"];


% 8. VRE
e2rmVRE = ["Quantico", "Woodbridge"; "Quantico", "Fredericksburg"; "Woodbridge", "Alexandria"; ...
    "Backlick Road_ 723", "Burke Centre_ 724"; "Alexandria", "Union Station, Washington, DC"];

% 9. Amtrak - ALL GOOD!

% 10. Baltimore MTA
e2rmBMR = ["Cherry Hill Light Rail Station_ 1518", "Mt Vernon Light Rail Station_ 1503"; "Camden Yards", "Penn Station, Baltimore"];

% 11. MBTA Subway
e2rmTSub = ["Braintree", "Quincy Center"; "Quincy Center", "JFK UMass"; "JFK UMass", "South Station"; "South Station", "Back Bay";...
    "Back Bay", "Ruggles"; "Forest Hills", "Ruggles"; "North Station", "Porter Square"; "North Station", "Malden Center"];

% 12. NYC Subway
e2rmNYC = ["Penn Station NYC", "Woodside"; "Grand Central Station", "125th St"; "Atlantic Ave/Barclays Center", "Nostrand Ave";...
    "Nostrand Ave", "Broadway"; "Broadway", "Jamaica"; "Forest Hills, NY", "Woodside"; "Woodside", "Hunters Point"];

% 13. PATH
e2rmPATH = ["Penn Station Newark", "Hoboken"];

% 14. SEPTA Metro
e2rmSEPM = ["30th St", "Pennsauken"; "15th St", "30th St"; "30th St", "15th St"];

% 15. DC Metro
e2rmDCM = ["Alexandria", "Franconia/Springfield"; "Alexandria", "Crystal City"; "Crystal City", "L'Enfants Plaza";...
    "Union Station, Washington, DC", "Alexandria"; "Union Station, Washington, DC", "L'Enfants Plaza"; "Union Station, Washington, DC", "New Carrollton";...
    "Silver Springs", "Union Station, Washington, DC"; "Union Station, Washington, DC", "Rockville"];

edgesToRm = {e2rmCTR, [], e2rmMARC, [], e2rmMNR, e2rmNJT, e2rmSEPR, e2rmVRE, [], e2rmBMR, e2rmTSub, e2rmNYC, e2rmPATH, e2rmSEPM, e2rmDCM};


%% 3. Build each subgraph

S = cell(size(pairs,1),1);

for k = 1:size(pairs,1)
    scale_k = pairs{k,1};
    syst_k  = pairs{k,2};
    include = false(n,1);
    for i = 1:n
        include(i) = any(nodeScales{i}  == scale_k) && ...
            any(nodeSystems{i} == syst_k);
    end
    Sk = subgraph(G, find(include));
    Sk.Edges.hubEnds = false(height(Sk.Edges),1);

    % Remove bad edges from other systems
    e2rmk = edgesToRm{k};
    if ~isempty(e2rmk)
        stk = reshape(findnode(Sk, e2rmk), size(e2rmk));
        Sk = rmedge(Sk, findedge(Sk, stk(:,1), stk(:,2)));
    end
    S{k} = Sk;
    pk = params;
    pk.numNodes = height(S{k}.Nodes);
    pk.G = S{k};
    if isfield(pk, "betweenness")
        pk.betweenness.NoN = pk.betweenness.NoN(include);
    end
    if isfield(pk, "rwBetweenness")
        pk.rwBetweenness.NoN = pk.rwBetweenness.NoN(include);
    end
    if isfield(pk, "degree")
        pk.degree.NoN =  pk.degree.NoN(include);
    end
    if isfield(pk, "closeness")
        pk.closeness.NoN =  pk.closeness.NoN(include);
    end
    if isfield(pk, "efficiencyDelta")
        pk.efficiencyDelta.NoN = pk.efficiencyDelta.NoN(include);
    end
    if k == 1
        paramsS = pk;
    else
        paramsS = [paramsS, pk];
    end
end

for k = 1:numel(S)
    Sk = S{k};
    for ki = 1:height(Sk.Edges)
        eki = Sk.Edges(ki,:);
        nski = findnode(Sk, eki.EndNodes);
        Tnski = Sk.Nodes(nski,:);

        if all(cellfun(@(x,y) (iscell(x) && ~isscalar(x)) || (iscell(y) && ~isscalar(y)), Tnski.system, Tnski.scale)) % if both EndNodes are interscale hubs connecting systems k, kj1, kj2...
            Sk.Edges.hubEnds(ki) = true;
        else
            Sk.Edges.hubEnds(ki) = false;
        end
    end
    S{k} = Sk;
end


end


%% Helper: convert value → column string vector, lower-case
function vec = toStringVec(val)
    if iscell(val)
        vec = lower(string(val(:)));   % cell --> column string array
    else
        vec = lower(string(val));      % char or string scalar
    end
end

