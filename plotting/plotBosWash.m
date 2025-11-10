function fig = plotBosWash(G, nodeCVals, params, fig, splt)

if nargin < 4 || isempty(fig)
    fig = figure;
else
    figure(fig)
end

if nargin > 4 && ~isempty(splt)
    subplot(splt)
end

% Get normalized geographic positions
xy_zscore = normalize([G.Nodes.lon, G.Nodes.lat]);

% Configure color map
nMapCols = 100;
if nargin >= 3
    cmap = params.graphPlot.colorMap(nMapCols);
else
    cmap = turbo(nMapCols);
end

if nargin < 2 || isempty(nodeCVals)
    dmgVals = zeros(numel(params.y0), 1);
    dmgVals(params.d_idx) = max(params.d,[],2);
    dmgValsMap = round((dmgVals - min(dmgVals))/(max(dmgVals)-min(dmgVals)) * (nMapCols-1) + 1);
    valsMap = dmgValsMap;
else
    if nargin > 2 && isfield(params.graphPlot, "cmapLimit") && ~isempty(params.graphPlot.cmapLimit)
        upperLim = params.graphPlot.cmapLimit(2);
        lowerLim = params.graphPlot.cmapLimit(1);
    else
        upperLim = max(nodeCVals);
        lowerLim = min(nodeCVals);
        valsMap = round((nodeCVals -lowerLim)/(upperLim-lowerLim) * (nMapCols-1) + 1);
    end
end

if all(isnan(valsMap))
    valsMap = ones(height(G.Nodes),1);
end

[valsSorted, vsi] = sort(valsMap, "descend");
colsSorted = cmap(valsSorted,:);
xySorted = xy_zscore(vsi,:);
GcScalePlot = plot(G, 'XData', xy_zscore(:,1), 'YData', xy_zscore(:,2), "NodeColor", cmap(valsMap,:), "EdgeColor", "k");

if isnumeric(params.graphPlot.markerSize)
    GcScalePlot.MarkerSize = params.graphPlot.markerSize;
elseif string(class(params.graphPlot.markerSize)) == "string" && params.graphPlot.markerSize == "proportional"
    mSize = rescale(nodeCVals, 2, 8);
    GcScalePlot.MarkerSize = mSize;
end


title("Northeast Rail Corridor", "Interpreter", "latex")
cb4 = colorbar;
cb4.Label.Interpreter = "latex";
cb4.TickLabelInterpreter = "latex";
cb4.Label.FontSize = 14;
if isfield(params.graphPlot, "colorBarLabelString")
    cb4.Label.String = params.graphPlot.colorBarLabelString;
end
colormap(fig, cmap)
ax5 = gca;
ax5.TickLabelInterpreter = "latex";
ax5.YTick = [];
ax5.XTick = [];
ax5.XAxis.Visible = "off";
ax5.YAxis.Visible = "off";

if nargin < 5
    fig.Position = [ 78 47 1220 1064];
else
    fig.Position = [ 30 271 1789 773];
end

pos = gca().Position;
xapf = @(x,pos,xl) pos(3)*(x-min(xl))/diff(xl)+pos(1); % 'x' Annotation Position Function
yapf = @(y,pos,yl) pos(4)*(y-min(yl))/diff(yl)+pos(2); % 'y' Annotation Position Function
annArr = cell(10,1);

switch params.graphPlot.nodeVar
    case "efficiency"
        annArr{1} = annotation('textarrow', xapf(xySorted(1,1)+[0.7,0], pos, xlim), yapf(xySorted(1,2)+[-1,-0.02], pos, ylim), 'String', ["\textbf{\#1: Newark Penn Station}"; "Newark, NJ"; "Regional, Intercity"], "Interpreter", "latex"); % "NJT, PATH, Amtrak";
        annArr{2} = annotation('textarrow', xapf(xySorted(2,1)+[-0.55,0],pos, xlim), yapf(xySorted(2,2)+[1,0.02], pos, ylim), 'String', ["\textbf{\#2: Liberty Int'l Airport}"; "Newark, NJ";  "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{3} = annotation('textarrow', xapf(xySorted(3,1)+[0.3,0.02],pos, xlim), yapf(xySorted(3,2)+[-0.7,-0.02], pos, ylim), 'String', ["\textbf{\#3: Metropark}"; "Iselin, NJ"; "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{4} = annotation('textarrow', xapf(xySorted(4,1)+[-0.55,0],pos, xlim), yapf(xySorted(4,2)+[0.65,0.0], pos, ylim), 'String', ["\textbf{\#4: New Brunswick}"; "New Brunswick, NJ";  "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{5} = annotation('textarrow', xapf(xySorted(5,1)+[-0.5,-0.03],pos, xlim), yapf(xySorted(5,2)+[0.39,0.02], pos, ylim), 'String', ["\textbf{\#5: Princeton Junction}"; "Princeton Junction, NJ";  "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{6} = annotation('textarrow', xapf(xySorted(6,1)+[0.1,0],pos, xlim), yapf(xySorted(6,2)+[-1,-0.1], pos, ylim), 'String', ["\textbf{\#6: Back Bay}"; "Boston, MA";  "Urban, Regional, Intercity"], "Interpreter", "latex"); % "MBTA, Amtrak";
        annArr{7} = annotation('textarrow', xapf(xySorted(7,1)+[0.7,0],pos, xlim), yapf(xySorted(7,2)+[-0.6,-0.0], pos, ylim), 'String', ["\textbf{\#7: Newark (Delaware)}"; "Newark, DE";  "Regional, Intercity"], "Interpreter", "latex"); % "MBTA, Amtrak";
        annArr{8} = annotation('textarrow', xapf(xySorted(8,1)+[0.5,0.03],pos, xlim), yapf(xySorted(8,2)+[-1,-0.05], pos, ylim), 'String', ["\textbf{\#8: Aberdeen}"; "Aberdeen, MD";  "Regional, Intercity"], "Interpreter", "latex"); % "MARC, Amtrak";
        annArr{9} = annotation('textarrow', xapf(xySorted(9,1)+[0.45,0.01],pos, xlim), yapf(xySorted(9,2)+[-1.45,-0.01], pos, ylim), 'String', ["\textbf{\#9: Baltimore Penn Station}"; "Baltimore, MD";  "Urban, Regional, Intercity"], "Interpreter", "latex"); % "BLR, MARC, Amtrak";
        annArr{10} = annotation('textarrow', xapf(xySorted(10,1)+[0.7,0.01],pos, xlim), yapf(xySorted(10,2)+[-0.5,-0.01], pos, ylim), 'String', ["\textbf{\#10: New York Penn Station}"; "Manhattan, NY"; "Urban, Regional, Intercity"], "Interpreter", "latex"); % "BLR, MARC, Amtrak";
        cb4.Label.String = "\textbf{Node importance: efficiency}";
    case "betweenness"
        annArr{1} = annotation('textarrow', xapf(xySorted(1,1)+[1,0], pos, xlim), yapf(xySorted(1,2)+[-.5,-0.01], pos, ylim), 'String', ["\textbf{\#1: New York Penn Station}"; "Manhattan, NY"; "Urban, Regional, Intercity"], "Interpreter", "latex"); % "NJT, PATH, Amtrak";
        annArr{2} = annotation('textarrow', xapf(xySorted(2,1)+[1.2,0], pos, xlim), yapf(xySorted(2,2)+[-1,-0.0], pos, ylim), 'String', ["\textbf{\#2: Newark Penn Station}"; "Newark, NJ"; "Regional, Intercity"], "Interpreter", "latex"); % "NJT, PATH, Amtrak";
        annArr{3} = annotation('textarrow', xapf(xySorted(3,1)+[0.34,0.02],pos, xlim), yapf(xySorted(3,2)+[-0.7,-0.02], pos, ylim), 'String', ["\textbf{\#3: Metropark}"; "Iselin, NJ"; "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{4} = annotation('textarrow', xapf(xySorted(4,1)+[-0.55,0],pos, xlim), yapf(xySorted(4,2)+[1.2,0.0], pos, ylim), 'String', ["\textbf{\#4: Liberty Int'l Airport}"; "Newark, NJ";  "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{5} = annotation('textarrow', xapf(xySorted(5,1)+[-0.55,0],pos, xlim), yapf(xySorted(5,2)+[0.85,0.0], pos, ylim), 'String', ["\textbf{\#6: New Brunswick}"; "New Brunswick, NJ";  "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{6} = annotation('textarrow', xapf(xySorted(6,1)+[-0.5,-0.02],pos, xlim), yapf(xySorted(6,2)+[0.5,0.02], pos, ylim), 'String', ["\textbf{\#7: Princeton Junction}"; "Princeton Junction, NJ";  "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{7} = annotation('textarrow', xapf(xySorted(7,1)+[-0.5,-0.03],pos, xlim), yapf(xySorted(7,2)+[0.3,0.02], pos, ylim), 'String', ["\textbf{\#5: Trenton Transit Center}"; "Trenton, NJ";  "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{8} = annotation('textarrow', xapf(xySorted(8,1)+[-0.1,-0.01],pos, xlim), yapf(xySorted(8,2)+[-0.7,-0.02], pos, ylim), 'String', ["\textbf{\#8: North Philadelphia}"; "Philadelphia, PA";  "Urban, Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{9} = annotation('textarrow', xapf(xySorted(9,1)+[0.6,-0.01],pos, xlim), yapf(xySorted(9,2)+[-0.8,-0.02], pos, ylim), 'String', ["\textbf{\#9: Cornwell Heights}"; "Cornwell Heights, PA";  "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        annArr{10} = annotation('textarrow', xapf(xySorted(10,1)+[-0.7,-0.01],pos, xlim), yapf(xySorted(10,2)+[1.6,0.01], pos, ylim), 'String', ["\textbf{\#10: New Rochelle}"; "New Rochelle, NY";  "Regional, Intercity"], "Interpreter", "latex"); % "NJT, Amtrak";
        cb4.Label.String = "\textbf{Node importance: betweenness}";
    case "none"
        % idk
end

if isfield(params.graphPlot, "textColor") && params.graphPlot.textColor
    for i = 1:numel(annArr)
        annArr{i}.TextColor = colsSorted(i,:);
    end
end

end