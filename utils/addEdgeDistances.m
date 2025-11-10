function G = addEdgeDistances(G, latVarName, lonVarName, distVarName)

% addEdgeDistances: compute Haversine length (km) of every edge in a graph

% G = ADD EDGEDISTANCES(G) takes an undirected or directed graph/digraph
% G whose node table contains columns:
%       lat  – latitude  in degrees
%       lon  – longitude in degrees

% It returns the same graph with a new/updated edge‑table column
% "Distance" that stores the great‑circle length (kilometres) of each
% edge, calculated with the vectorised Haversine formula.

% The routine is fully vectorised—no loops—so it scales to millions of
% edges. Runtime is dominated by the trigonometric calls (~5 x |E| FLOP).

% ------------------------------------------------------------------------
% 1. Fill any missing inputs with defaults
% ------------------------------------------------------------------------
% If the variable names for node lat, lon are not given, assume they are
% "lat" and "lon". If a name for the edge distance variable is not given,
% name the variable "havdist_km"
if nargin == 1
    latVarName = "lat";
    lonVarName = "lon";
    distVarName = "havdist_km";
end

% ------------------------------------------------------------------------
% 2. Resolve edge endpoints to numeric node indices
% ------------------------------------------------------------------------
E = G.Edges.EndNodes; % |E| x 2 node‑index array

if isnumeric(E) % already indices
    s = E(:,1);
    t = E(:,2);
else % cellstr/string
    % "findnode" is vectorised and compiled (cheap)
    s = findnode(G, E(:,1));
    t = findnode(G, E(:,2));
end

% ------------------------------------------------------------------------
% 3. Vectorised Haversine on the edge endpoints
% ------------------------------------------------------------------------
R      = 6371.0088; % mean Earth radius, kilometers
latRad = deg2rad(G.Nodes.(latVarName)(:)); % m x 1, in radians
lonRad = deg2rad(G.Nodes.(lonVarName)(:));

lat1 = latRad(s);
lat2 = latRad(t);
lon1 = lonRad(s);
lon2 = lonRad(t);

dLat = lat2 - lat1;
dLon = lon2 - lon1;

a = sin(dLat/2).^2 + cos(lat1).*cos(lat2).*sin(dLon/2).^2;
c = 2 .* atan2(sqrt(a), sqrt(1-a));

G.Edges.(distVarName) = R * c; % |E| x 1 numeric column (km)

end