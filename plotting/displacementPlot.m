function plt = displacementPlot(x1, y1, x2, y2, cols, ttl, axLbls, mkrShp)

if nargin < 8 || isempty(mkrShp)
    mkrShp = "o";
end

% Plot lines
for ki = 1:numel(x1)
    plot([x1(ki), x2(ki)]', [y1(ki), y2(ki)]', LineWidth=5.3, Color=[cols(ki,:) 0.2])
end

scatter(x2, y2, 40, cols, mkrShp, MarkerFaceAlpha=0.8, LineWidth=1);
scatter(x1, y1, 50, cols, mkrShp, "filled", MarkerFaceAlpha=0.8);

if nargin > 5 && ~isempty(ttl)
    title(ttl, Interpreter="latex", FontSize=14)
end

if nargin > 6 && ~isempty(axLbls)
    if ~isempty(axLbls(1)) && axLbls(1) ~= "" % plot x axis label if included
        xlabel(axLbls(1), Interpreter="latex", FontSize=14)
    end
    if ~isempty(axLbls(2)) && axLbls(2) ~= ""
        ylabel(axLbls(2), Interpreter="latex", FontSize=14)
    end
end

plt = gca;

end