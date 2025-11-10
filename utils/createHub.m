function Gprime = createHub(G, nodeNames, hubName, hubMethod)

% createHub.m creates a hub node by replacing or joining specified nodes, with geographic midpoint

% Gprime = CREATEHUB(G, nodeNames, hubName, hubMethod) takes an undirected
% graph G with node attributes Name, lat, lon, and index. nodeNames can be:
%     - cell array of node names (char or string)
%     - numeric array matching G.Nodes.index values
% hubName is the new node's name (string or char), and hubMethod is
% either 'replace' or 'join'. The new hub node is assigned the geographic
% midpoint (mean latitude and longitude) of the specified nodes

% 'replace': remove nodeNames and connect hub to all their unique neighbors
% 'join': remove edges among nodeNames and connect hub to each

% Combination of hand-coding (JRW) and ChatGPT-5 Pro

% Validate that hubName does not already exist
if any(strcmp(G.Nodes.Name, hubName))
    error('createHub:HubExists', 'The hub name "%s" already exists.', hubName);
end


% Determine node indices from nodeNames input
if isnumeric(nodeNames)
    % Interpret numeric entries as values of G.Nodes.index
    idxNodes = arrayfun(@(x) find(G.Nodes.index==x, 1), nodeNames, "UniformOutput", false);
    if any(cellfun(@isempty, num2cell(idxNodes)))
        missing = nodeNames(cellfun(@isempty, num2cell(idxNodes)));
        error('createHub:IndexNotFound', 'No node with index value(s): %s', strjoin(string(missing), ', '));
    end
    if iscell(idxNodes)
        idxNodes = cell2mat(idxNodes);
    end
    nodeNames = G.Nodes.Name(idxNodes);
else
    if isstring(nodeNames) % Convert string/char inputs to cellstr
        nodeNames = cellstr(nodeNames);
    elseif ischar(nodeNames)
        nodeNames = {nodeNames};
    end
    [found, idxNodes] = ismember(nodeNames, G.Nodes.Name); % Map names to indices
    if ~all(found)
        missing = nodeNames(~found);
        error('createHub:NodeNotFound', 'Nodes not found: %s', strjoin(missing, ', '));
    end
end

latMid = mean(G.Nodes.lat(idxNodes)); % Compute geographic midpoint
lonMid = mean(G.Nodes.lon(idxNodes));

hubMethod = char(hubMethod); % Validate hubMethod
if ~ismember(hubMethod, {'replace','join'})
    error('createHub:InvalidMethod', 'hubMethod must be "replace" or "join".');
end

switch hubMethod
    case 'replace'
        % Collect all unique neighbors of specified nodes
        A = adjacency(G);
        neighborMask = any(A(idxNodes,:), 1);
        nbrIdx = find(neighborMask);
        nbrIdx = setdiff(nbrIdx, idxNodes);

        % Add the new hub node
        G1 = addnode(G, hubName);
        hubIdx = findnode(G1, hubName);
        % Assign midpoint coordinates
        G1.Nodes.lat(hubIdx) = latMid;
        G1.Nodes.lon(hubIdx) = lonMid;
        
        % Handle system and scale variables
        sysExpanded = G.Nodes.system(idxNodes);
        sysExpanded = vertcat(sysExpanded{:});
        if ischar(sysExpanded)
            sysExpanded = cellstr(sysExpanded);
        end
        sysExpanded(cellfun(@(x) isempty(x), sysExpanded)) = [];
        hubSystems = unique(sysExpanded);
        scaleExpanded = G.Nodes.scale(idxNodes);
        scaleExpanded = vertcat(scaleExpanded{:});
        if ischar(scaleExpanded)
            scaleExpanded = cellstr(scaleExpanded);
        end
        scaleExpanded(cellfun(@(x) isempty(x), scaleExpanded)) = [];
        hubScales  = unique(scaleExpanded);
        G1.Nodes.system(hubIdx) = {hubSystems};
        G1.Nodes.scale(hubIdx) = {hubScales};
        G1.Nodes.rshp(hubIdx) = sum(G.Nodes.rshp(idxNodes));
       

        % Remove the original nodes
        G1 = rmnode(G1, idxNodes);

        % Connect hub to each unique neighbor by name
        neighborNames = G.Nodes.Name(nbrIdx);
        Gprime = addedge(G1, repmat({hubName}, numel(neighborNames),1), neighborNames, ones(1,numel(neighborNames)));

    case 'join'
        % Remove all edges among the specified nodes
        G1 = G;
        pairs = nchoosek(idxNodes,2);
        for k = 1:size(pairs,1)
            if findedge(G1, pairs(k,1), pairs(k,2)) > 0
                G1 = rmedge(G1, pairs(k,1), pairs(k,2));
            end
        end

        % Add the new hub node
        G1 = addnode(G1, hubName);
        hubIdx = findnode(G1, hubName);
        % Assign midpoint coordinates
        G1.Nodes.lat(hubIdx) = latMid;
        G1.Nodes.lon(hubIdx) = lonMid;

        % Connect hub to each of the original nodes by index
        Gprime = addedge(G1, repmat(hubIdx, numel(idxNodes),1), idxNodes);
end
end
