function G = initBosWashGraph(fpath)

if nargin >= 1 && ~isempty(fpath)
    BWG = load(fpath);
else
    BWG = load("data\BosWash_graph_26Sep2025.mat"); % default location
end

G = BWG.G;

% Add edge weights as haversine distances
G = addEdgeDistances(G, "lat", "lon", "Weight");
G.Edges.Weight = rescale(G.Edges.Weight, 1, 2);

end