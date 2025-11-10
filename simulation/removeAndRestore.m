function Gstates = removeAndRestore(G, B)

% removeAndRestore  Remove then re-add nodes by name, capturing each graph state,
%                   and preserving edge weights (if present).
%   Gstates = removeAndRestore(G, B):

% IN: 
% G — a graph or digraph whose Nodes table has a 'Name' variable
% B — an N×1 cell array of node names to remove/add in order

% OUT: 
% 1×(2N+1) cell array of graph objects:
%   Gstates{1} = original G
%   Gstates{k+1} for k=1:N = after removing B{k}
%   Gstates{N+1+k} for k=1:N = after adding back B{k}, with original edges & weights

% Store original graph & full edge table
G0 = G;
origE   = G0.Edges; % table with EndNodes (and maybe Weight, etc.)
origEN  = origE.EndNodes; % E×2 cell array (string or cellstr)
if isstring(origEN)
    origEN = cellstr(origEN);
end
hasW    = ismember('Weight', origE.Properties.VariableNames);

N = numel(B);
Gstates = cell(1, 2*N + 1);
step = 1;
Gstates{step} = G; % step 0: original

% 1. remove each node in sequence -----------------------------------------
for k = 1:N
    name = B{k};
    idx  = findnode(G, name);
    if idx==0
        warning('Node "%s" not found; skipping removal.', name);
    else
        G = rmnode(G, idx);
    end
    step = step + 1;
    Gstates{step} = G;
end

% 2. add back nodes & restore edges (+weights) ----------------------------
for k = 1:N
    name = B{k};
    G = addnode(G, {name});

    % restore any original edge (name <--> other) if "other" still exists
    rows1 = find(strcmp(origEN(:,1), name));
    for r = rows1'
        other = origEN{r,2};
        if findnode(G, other)~=0
            if hasW
                w = origE.Weight(r);
                G = addedge(G, name, other, w);
            else
                G = addedge(G, name, other);
            end
        end
    end

    rows2 = find(strcmp(origEN(:,2), name));
    for r = rows2'
        other = origEN{r,1};
        if findnode(G, other)~=0
            if hasW
                w = origE.Weight(r);
                G = addedge(G, name, other, w);
            else
                G = addedge(G, name, other);
            end
        end
    end

    step = step + 1;
    Gstates{step} = G;
end
end