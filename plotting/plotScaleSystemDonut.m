function plotScaleSystemDonut(G, fig, splt)
% plotScaleSystemDonut plots a nested donut chart of node counts by scale and system.
% Uses LaTeX multi-line labels directly on slices and a modern interpolated color palette.

% Combination of hand-coding (JRW) and ChatGPT 5 Pro

% Extract scale categories
for k = 1:height(G.Nodes)
    sck = G.Nodes.scale{k};
    if ischar(sck) || (isstring(sck) && isscalar(sck))
        G.Nodes.Scale{k} = strrep(sck, 'National', 'Intercity');
    elseif iscell(sck)
        G.Nodes.Scale{k} = cellfun(@(x) strrep(x, 'National', 'Intercity'), sck, UniformOutput=false);
    end
end

scaleCats = {'Urban','Commuter','National'};
scCombs = nchoosek(scaleCats, 2);
combPairs = cellfun(@(x,y) {x,y}, scCombs(:,1), scCombs(:,2), "UniformOutput", false);
scaleCats = [scaleCats, combPairs', {scaleCats}];
scaleStrs = cellfun(@(x) sort(string(x), 'descend'), G.Nodes.scale, 'UniformOutput', false);

numScale   = numel(scaleCats);
countsScale = zeros(1,numScale);
systemCats  = cell(1,numScale);
countsSystem= cell(1,numScale);

for i = 1:numScale
    if iscell(scaleCats{i})
        idx = false(height(G.Nodes),1);
        sci = string(scaleCats{i});
        for j = 1:numel(scaleStrs)
            if numel(sci) == numel(scaleStrs{j})
                idx(j) = all(strcmp(sort(sci', 'descend'), scaleStrs{j}));
            else
                idx(j) = false;
            end
        end
    else
        idx = strcmp(G.Nodes.scale, scaleCats{i});
    end
    countsScale(i) = sum(idx);
    if iscell(scaleCats{i})
        if numel(sci) == 2
            systemCats{i} = {strjoin(scaleCats{i}, ' & ')};
        elseif numel(sci) == 3
            systemCats{i} = {'Urban, Regional, & Intercity'};
        end
    else
        systemCats{i}  = unique(G.Nodes.system(idx));
    end
    countsSystem{i} = zeros(1,numel(systemCats{i}));
    if isscalar(systemCats{i}) 
        countsSystem{i}(1) = sum(idx);
    else
        for j = 1:numel(systemCats{i})
            countsSystem{i}(j) = sum(idx & strcmp(G.Nodes.system, systemCats{i}{j}));
        end
    end
end

% Define palettes (dark to light) for each scale
palettes = { ...
    [0.6,0,0],   [1,0.6,0.6]; ...  % Urban: dark red --> light red
    [0.4,0,0.6], [0.8,0.6,1]; ...  % Commuter: dark purple --> light purple
    [0,0,0.6],   [0.6,0.6,1]; ...       % National: dark blue --> light blue
    hex2rgb("#f1bb66"), hex2rgb("#f1bb66"); ... % Urban & Commuter
    hex2rgb("#efb04d"), hex2rgb("#efb04d");... % Urban & National
    hex2rgb("#e79f2e") , hex2rgb("#e79f2e"); ... % Commuter & National
    hex2rgb("#d08c22"), hex2rgb("#d08c22") % All 3 scales
};
innerColors = zeros(numScale,3);
outerColors = cell(1,numScale);
for i = 1:numScale
    cDark = palettes{i,1};
    cLight= palettes{i,2};
    innerColors(i,:) = (cDark + cLight)/2;
    m = numel(systemCats{i});
    outerColors{i} = [ ...
        linspace(cDark(1),cLight(1),m)', ...
        linspace(cDark(2),cLight(2),m)', ...
        linspace(cDark(3),cLight(3),m)'  ...
    ];
end

% Create figure
if nargin < 2 || isempty(fig)
    figure;
else
    figure(fig)
end

if nargin > 2 && ~isempty(splt)
    subplot(splt)
end

axis equal off;
hold on;

% Radii definitions
rInner0 = 0;    rInner1 = 0.3;
rOuter0 = 0.35; rOuter1 = 1;

totalNodes = sum(countsScale);
thetaInner = [0, cumsum(countsScale/totalNodes)*2*pi];
scaleCatLbls = {'Urban', 'Regional', 'Intercity', 'Urban & Regional', 'Urban & Intercity', 'Regional & Intercity', 'Urban, Regional, & Intercity'};
hTxtArr = cell(numScale,1);
% Plot inner ring
for i = 1:numScale
    t = linspace(thetaInner(i), thetaInner(i+1), 200);
    % draw slice
    patch([rInner0*cos(t), fliplr(rInner1*cos(t))], ...
        [rInner0*sin(t), fliplr(rInner1*sin(t))], ...
        innerColors(i,:), 'EdgeColor','none');
    % label slice
    mid     = (thetaInner(i)+thetaInner(i+1))/2;
    xm      = (rInner0+rInner1)/2 * cos(mid);
    ym      = (rInner0+rInner1)/2 * sin(mid);
    pct     = countsScale(i)/totalNodes*100;
    lbl     = ["$\textbf{" + string(scaleCatLbls{i}) + "}$"; ...
        "$\mathbf{" + num2str(countsScale(i)) + " \textrm{ } (" + num2str(round(pct,1)) + "\%)}$"];
    lbl(1) = strrep(lbl(1), "&", "\&");
    if  string(scaleCats{i}) == "National"
        txtCol = "w";
    else
        txtCol = "k";
    end
    hTxt    = text(xm, ym, lbl, 'HorizontalAlignment','center', ...
        'Interpreter','latex', "Color", txtCol, 'FontSize', 16);
    hTxtArr{i} = hTxt;
    % ensure text is on top
    %uistack(hTxt, 'bottom');
    uistack(hTxt, 'top');
end

for ii = 1:numel(hTxtArr)
    uistack(hTxtArr{ii}, 'top')
end


hTxtArr = cell(numScale*sum(cellfun(@(x) numel(x), countsSystem)),1);
hTxtIdx = 1;

% Plot outer ring
systemCatLbls = cell(1,numel(systemCats));
systemCatLbls{1} = {'MTA (Baltimore)'; 'MBTA (Boston)'; 'MTA (New York)'; ...
    'PATH (New York/New Jersey)'; 'SEPTA (Philadelphia)'; 'WMATA (Washington, DC)'};
systemCatLbls{2} = {'CT Rail (Connecticut)'; 'Long Island Rail Road (New York)'; 'Maryland Area Rail Commuter';...
    'MBTA Commuter Rail (Massachusetts)'; 'Metro-North (New York)'; 'New Jersey Transit';...
    'SEPTA (Philadelphia)'; 'Virginia Railway Express'};
systemCatLbls{3} = {'Amtrak'};
systemCatLbls(4:end) = {{'Urban & Regional'}, {'Urban & Intercity'}, {'Regional & Intercity'}, {'Urban, Regional, & Intercity'}};


for i = 1:numScale
    cumFrac = [0, cumsum(countsSystem{i}/countsScale(i))];
    for j = 1:numel(countsSystem{i})
        t0 = thetaInner(i) + cumFrac(j)*(thetaInner(i+1)-thetaInner(i));
        t1 = thetaInner(i) + cumFrac(j+1)*(thetaInner(i+1)-thetaInner(i));
        t  = linspace(t0, t1, 200);
        % draw slice
        patch([rOuter0*cos(t), fliplr(rOuter1*cos(t))], ...
              [rOuter0*sin(t), fliplr(rOuter1*sin(t))], ...
              outerColors{i}(j,:), 'EdgeColor','none');
        % label slice
        mid  = (t0 + t1)/2;
        xm   = (rOuter0+rOuter1)/2 * cos(mid);
        ym   = (rOuter0+rOuter1)/2 * sin(mid);
        pct  = countsSystem{i}(j)/totalNodes*100;
        lbl  = ["$\textbf{" + string(systemCatLbls{i}{j}) + "}$"; ...
                "$\mathbf{" + num2str(countsSystem{i}(j)) + " \textrm{ } (" + num2str(round(pct,1)) + "\%)}$"];
        lbl(1) = strrep(lbl(1), "&", "\&");
        if any(string(systemCats{i}{j}) == ["LIRR", "MARC", "MBTA"]) && string(scaleCats{i}) == "Commuter"
            txtCol = "w";
        elseif string(systemCats{i}{j}) == "Baltimore" && string(scaleCats{i}) == "Urban"
            txtCol = "w";
        else
            txtCol = "k";
        end
         hTxt = text(xm, ym, lbl, 'HorizontalAlignment','center', ...
                     'Interpreter','latex', "Color", txtCol);
         hTxtArr{hTxtIdx} = hTxt;
        hTxtIdx = hTxtIdx + 1;
        uistack(hTxt, 'top');
    end
end

for ii = 1:numel(hTxtArr)
    uistack(hTxtArr{ii}, 'top')
end

% Title
title('$\textbf{Node counts by scale and system}$', ...
      'Interpreter','latex', 'FontSize',14);

hold off;

end