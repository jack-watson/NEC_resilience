
%% Main (top-level) script: megaregional RTN analysis
% Title: "Megaregional adaptation andrecovery priorities preserve local 
% resilience in the US Northeast rail corridor"
% Jack R. Watson, Samrat Chatterjee, and Auroop R. Ganguly, 2025

% All code produced and owned by Jack R. Watson
% Contact: watson.jac@northeastern.edu / jwatson294@gmail.com
% GitHub: https://github.com/jack-watson

addpath(genpath(cd)) % add subfolders to path
bw_fpath = "data\BosWash_graph_26Sep2025.mat"; % specify network data path
G = initBosWashGraph(bw_fpath, false); % load NEC graph object

% [fList,pList] = matlab.codetools.requiredFilesAndProducts('bosWashMain.m');
% fList = fList';

%% Configure parameters

% set this to TRUE to produce miscellanious visualizations used for
% mid-analysis diagnostics, descriptive stats, etc.
% WARNING: turning this on will make this script output a LOT of plots
miscPlotsTF = true;

params.G = G; % store graph info in params structure
params.numNodes = height(G.Nodes);
params.wgtFun = @(Gx) rescale(centrality(Gx, "betweenness"), 0.01, 0.99);
params.isWeighted = true; % true if edges are distance-weighted, false if not
params.sysPerformanceFun = @(y, pms) networkEfficiency(y, pms); % function that measures performance

% For NEC-specific plotting later on
xy_zscore = normalize([G.Nodes.lon, G.Nodes.lat]);
scale = G.Nodes.scale;
urbanTF = strcmp(scale, 'Urban');
commTF  = strcmp(scale, 'Commuter');
natTF   = strcmp(scale, 'National');
scaleCols = zeros(length(scale),3);
scaleCols(urbanTF,:) = ones(sum(urbanTF),3).*[0.85,0.33,0.10];
scaleCols(commTF,:)  = ones(sum(commTF),3).*[0.72,0.27,1.00];
scaleCols(natTF,:)   = ones(sum(natTF),3).*[0.00,0.45,0.74];

%% Assign ridership/OD with gravity model, calculate node importance

% Expect gravityMdl.m to take several minutes to run
[gmT, gmEF, gmInfo] = gravityMdl(G, true, true);
% Calculate ridership-weighted betweenness centrality using gravMdl output
[rwBC, info] = rwbc(G, gmT, false); % Also computationally intensive
G.Nodes.rwbc = rwBC; % Add rwbc as a node attribute in G

% Get global network efficiency upon removing nodes one by one
% Can take a few minutes to run on a consumer-grade workstation
params.resOpts.weightByRshpTF = false;
nodeRmEff = calcNodeRmImportance(G, params, @networkEfficiency);
nodeRmEff = rescale(1 - nodeRmEff, 0, 1);

% Calculate betweenness centrality
bc = centrality(params.G, 'betweenness', Cost = G.Edges.Weight);
[~, bcInd] = sort(bc, "descend");

%% 2.1 Descriptive visualization of BosWash corridor network

params.graphPlot.colorMap = @turbo;
params.graphPlot.markerSize = "proportional";

if miscPlotsTF
    bwFigTest = figure("Name", "Boston-Washington viz: basic");
    params.graphPlot.nodeVar = "none";
    params.graphPlot.textColor = false;
    plotBosWash(G, nodeRmEff, params, bwFigTest)

    bwFig1 = figure("Name", "Boston-Washington viz: efficiency");
    params.graphPlot.nodeVar = "efficiency";
    params.graphPlot.textColor = false;
    plotBosWash(G, nodeRmEff, params, bwFig1)

    bwFig2 = figure("Name", "Boston-Washington viz: betweenness");
    params.graphPlot.nodeVar = "betweenness";
    params.graphPlot.textColor = false;
    plotBosWash(G, bc, params, bwFig2)

    bwFig3 = figure("Name", "Boston-Washington viz: donut");
    plotScaleSystemDonut(G, bwFig3)

    bwFigR = figure("Name", "Boston-Washington viz: ridership");
    params.graphPlot.nodeVar = "none";
    params.graphPlot.textColor = false;
    plotBosWash(G, G.Nodes.rshp, params, bwFigR)

    bwFigR2 = figure("Name", "Boston-Washington viz: ridership-weighted betweenness centrality");
    params.graphPlot.nodeVar = "none";
    params.graphPlot.textColor = false;
    plotBosWash(G, G.Nodes.rwbc, params, bwFigR2)
end

%% Measure node-level attributes in isolation (local) and in NoN (global)

% Calculate global measures of importance (S-G)
params.betweenness.NoN = bc;
params.rwBetweenness.NoN = rwBC;
params.resOpts.weightByRshpTF = false;
params.efficiencyDelta.NoN = 1 - (calcNodeRmImportance(G, params, @networkEfficiency) ./ networkEfficiency(ones(height(params.G.Nodes),1), params));
params.closeness.NoN = centrality(params.G, 'closeness', Cost = G.Edges.Weight);
params.degree.NoN = centrality(params.G, 'degree');

% Disaggregate individual RTNs for local analysis (S-L)
[S, pS] = splitByScaleSys(G, params);
rshpVarNames = "rshp_" + ["CTRail", "LIRR", "MARC", "MBTACom", "MNR", "NJT", ...
    "SEPTAReg", "VRE", "Amtrak", "Baltimore", "T", "NYC", "PATH", "SEPTAMetro", "WMATA"];

% Graph plots of each individual RTN
if miscPlotsTF
    for i = 1:numel(S)
        pS(i).graphPlot.nodeVar = rshpVarNames(i);
        bwFig = plotBosWash(S{i}, ones(height(S{i}.Nodes),1), pS(i));
        bwFig.Children(2).Children.NodeLabel = [];
    end
end

for k = 1:numel(S) % iterate through each RTN in the NEC...
    
    gk = S{k}; % gk is subgraph (RTN) k

    % Get OD matrix from gravity model
    rshpk = gk.Nodes.(rshpVarNames(k));
    rshpk(isnan(rshpk)) = 1;
    [ODk, EFk, gmInfok] = gravityMdl(gk, true, true, rshpk);

    % Measure local (S-L) importance: DC, CC, DE, BC, RWBC
    bck = centrality(gk, 'betweenness', Cost = gk.Edges.Weight);
    [rwbck, ~] = rwbc(gk, ODk, false);
    cck = centrality(gk, 'closeness', Cost = gk.Edges.Weight);
    edk = 1 - (calcNodeRmImportance(gk, pS(k), @networkEfficiency) ./ networkEfficiency(ones(height(gk.Nodes),1), pS(k)));
    dck = centrality(gk, 'degree');

    pS(k).betweenness.indv = bck; % store raw measurements in params structs
    pS(k).closeness.indv = cck;
    pS(k).efficiencyDelta.indv = edk;
    pS(k).degree.indv = dck;
    pS(k).rwBetweenness.indv = rwbck;
    pS(k).ODMatrix.indv = ODk;

    bckNorm = rescale(bck, 0, 1); % normalize all measures to [0,1]
    rwbckNorm = rescale(rwbck, 0, 1);
    bcNoNkNorm = rescale(pS(k).betweenness.NoN, 0, 1);
    rwbcNoNkNorm = rescale(pS(k).rwBetweenness.NoN, 0, 1);
    cckNorm = rescale(cck, 0, 1);
    ccNoNkNorm = rescale(pS(k).closeness.NoN, 0, 1);

    pS(k).betweenness.ratio =  bckNorm./bcNoNkNorm; % ratio of measures for S-L/S-G
    pS(k).rwBetweenness.ratio = rwbckNorm./rwbcNoNkNorm;
    pS(k).closeness.ratio = cckNorm./ccNoNkNorm;
    pS(k).degree.ratio = dck./pS(k).degree.NoN;
    pS(k).efficiencyDelta.ratio = edk./pS(k).efficiencyDelta.NoN;

    for i = 1:numel(bck) % if measure is 0 for both S-L and S-G, set to 1 to avoid Inf or NaN
        if bckNorm(i) == 0 && bcNoNkNorm(i) == 0
            pS(k).betweenness.ratio(i) = 1;
        end
        if rwbckNorm(i) == 0 && rwbcNoNkNorm(i) == 0
            pS(k).rwBetweenness.ratio(i) = 1;
        end
        if cckNorm(i) == 0 && ccNoNkNorm(i) == 0
            pS(k).closeness.ratio(i) = 1;
        end
        if edk(i) == 0 && pS(k).efficiencyDelta.NoN(i) == 0
            pS(k).efficiencyDelta.ratio(i) = 1;
        end
        if dck(i) == 0 && pS(k).degree.NoN(i) == 0
            pS(k).degree.ratio(i) = 1;
        end
    end

    % store results as node attributes for each RTN graph Sk
    [S{k}.Nodes.btw_ratio, pS(k).G.Nodes.btw_ratio] = deal(pS(k).betweenness.ratio);
    [S{k}.Nodes.rwbtw_ratio, pS(k).G.Nodes.rwbtw_ratio] = deal(pS(k).rwBetweenness.ratio);
    [S{k}.Nodes.deg_ratio, pS(k).G.Nodes.deg_ratio] = deal(pS(k).degree.ratio);
    [S{k}.Nodes.cls_ratio, pS(k).G.Nodes.cls_ratio] = deal(pS(k).closeness.ratio);
    [S{k}.Nodes.eff_ratio, pS(k).G.Nodes.eff_ratio] = deal(pS(k).efficiencyDelta.ratio);

end

%% Diagnostic plots of importance measures

pIdx = 4; % set to select which Sk to visualize (in range [1,15])

if miscPlotsTF
    bwBtwRatioFig = figure("Name", "Boston-Washington viz: ratio of indv to NoN betweenness");
    params.graphPlot.nodeVar = "none";
    params.graphPlot.textColor = false;
    logBtwRatio = log(pS(pIdx).betweenness.ratio);
    logBtwRatio(logBtwRatio == -Inf) = 0;
    plotBosWash(S{pIdx}, pS(pIdx).betweenness.ratio, pS(pIdx), bwBtwRatioFig)
    clim([min(pS(pIdx).betweenness.ratio), max(pS(pIdx).betweenness.ratio)])

    bwDegRatioFig = figure("Name", "Boston-Washington viz: ratio of indv to NoN degree");
    degRatio = pS(pIdx).degree.ratio;
    degRatio(degRatio == -Inf) = 0;
    plotBosWash(S{pIdx}, degRatio, pS(pIdx), bwDegRatioFig)
    clim([min(pS(pIdx).degree.ratio), max(pS(pIdx).degree.ratio)])
end

%% Visualize change in prioritization sequence

tK = 10; % number of top-K ranked stations to visualize; set (tK = 10, tKallTF = false) to replicate Figure 3 & Figure 4
tKallTF = true; % overrides tK, sets tK = |Vk| for subgraph Sk(Vk,Ek). Set to TRUE to replicate Figure 5 & Figure 6
alluvTF = true; % enable/disable alluvial plots; set to TRUE to replicate Figure 4 (WARNING: will output many plots)

dr = zeros(numel(S),1);
fig1 = figure; tl1 = tiledlayout(fig1, 4,4);
fig2 = figure; tl2 = tiledlayout(fig2, 4,4);
fig3 = figure; tl3 = tiledlayout(fig3, 3, 1);
fig4 = figure; tl4 = tiledlayout(fig4, 4,4);
urbColsQuad = [hex2rgb("#f5b7b1"); hex2rgb("#d98880"); hex2rgb("#e74c3c"); hex2rgb("#b03a2e"); hex2rgb("#822b21")]; 
regColsQuad = [hex2rgb("#d5f5e3"); hex2rgb("#82e0aa"); hex2rgb("#52be80"); hex2rgb("#1e8449"); hex2rgb("#166034")];
intColsQuad = [hex2rgb("#d6eaf8"); hex2rgb("#7fb3d5"); hex2rgb("#2e86c1"); hex2rgb("#1a5276"); hex2rgb("#0d293b")];
quadColsIndv = [hex2rgb("#273746"); hex2rgb("#82e0aa"); hex2rgb("#7fb3d5"); hex2rgb("#8e44ad")];
quadColsNoN = [hex2rgb("#f4d03f"); hex2rgb("#eb984e"); hex2rgb("#e74c3c"); hex2rgb("#db36a7")];

sysNames = ["CT Rail", "LIRR", "MARC", "MBTA Commuter", "Metro-North", ...
    "NJ Transit", "SEPTA Regional", "VRE", "Amtrak", "MTA (Baltimore)", ...
    "MBTA Subway", "NYC Subway", "PATH", "SEPTA Metro", "Washington Metro"];
scaleCols = [hex2rgb("#cb4335 "); hex2rgb("#27ae60"); hex2rgb("#007cff")]; % red = urban, green = regional, blue = intercity
Ccol = hex2rgb("#FFC300");
tileOrder = [1 2 5 6 ...
            9 10 13 14 ...
            15 3 4 7 ...
            8 11 12 16];
alluvOpts = struct; alluvOpts.textLabels = false; alluvOpts.rankFontSize = 13;
[E_arr_bc, E_arr_rwbc, E_arr_dc, E_arr_cc, E_arr_ed, ...
    dChg_arr_bc, dChg_arr_rwbc, dChg_arr_dc, dChg_arr_cc, dChg_arr_ed] = deal(cell(numel(S), 2));

textSz = 14; % fontSize arg in xlabel/ylabel/legend/title calls 
patchColors = cell(3,1); pcLen = 15;
pcUstart = hex2rgb("#C90000"); pcUend = hex2rgb("#F2D3D3");
pcRstart = hex2rgb("#0BA100"); pcRend = hex2rgb("#D4F2D3");
pcIstart = hex2rgb("#000DC9"); pcIend = hex2rgb("#D3D5F2");
patchColors{1} = [linspace(pcUstart(1),pcUend(1),pcLen)', linspace(pcUstart(2),pcUend(2),pcLen)', linspace(pcUstart(3),pcUend(3),pcLen)'];
patchColors{2} = [linspace(pcRstart(1),pcRend(1),pcLen)', linspace(pcRstart(2),pcRend(2),pcLen)', linspace(pcRstart(3),pcRend(3),pcLen)'];
patchColors{3} = [linspace(pcIstart(1),pcIend(1),pcLen)', linspace(pcIstart(2),pcIend(2),pcLen)', linspace(pcIstart(3),pcIend(3),pcLen)'];
resOpts = struct; resOpts.weightByRshpTF = false;

for i = 1:numel(S) 
    
    % Get importance measures for S-L and S-G =============================
    [bcNoN, bcNoNInd]       = sort(pS(i).betweenness.NoN, "descend");
    [bcIndv, bcIndvInd]     = sort(pS(i).betweenness.indv, "descend");
    [rwbcNoN, rwbcNoNInd]   = sort(pS(i).rwBetweenness.NoN, "descend");
    [rwbcIndv, rwbcIndvInd] = sort(pS(i).rwBetweenness.indv, "descend");
    [efNoN, efNoNInd]       = sort(pS(i).efficiencyDelta.NoN, "descend");
    [efIndv, efIndvInd]     = sort(pS(i).efficiencyDelta.indv, "descend");
    [dgNoN, dgNoNInd]       = sort(pS(i).degree.NoN, "descend");
    [dgIndv, dgIndvInd]     = sort(pS(i).degree.indv, "descend");
    [clNoN, clNoNInd]       = sort(pS(i).closeness.NoN, "descend");
    [clIndv, clIndvInd]     = sort(pS(i).closeness.indv, "descend");
    resOpts.ODMatrix = pS(i).ODMatrix;
    resOpts.ODMatrix.NoN = gmT;
    resOpts.ODMatrix
    resOpts.ODMatrix.sum.NoN = sum(gmT, "all");
    resOpts.ODMatrix.sum.indv = sum(resOpts.ODMatrix.indv, "all");
    if tKallTF
        tK = numel(bcNoN);
    end


    % Get station names for top K-ranked nodes ============================
    bcTop10Ind  = [bcIndvInd(1:tK), bcNoNInd(1:tK)];
    bcNoNNames  = pS(i).G.Nodes.Name(bcTop10Ind(:,2));
    bcIndvNames = pS(i).G.Nodes.Name(bcTop10Ind(:,1));
    rwbcTop10Ind  = [rwbcIndvInd(1:tK), rwbcNoNInd(1:tK)];
    rwbcNoNNames  = pS(i).G.Nodes.Name(rwbcTop10Ind(:,2));
    rwbcIndvNames = pS(i).G.Nodes.Name(rwbcTop10Ind(:,1));
    clTop10Ind  = [clIndvInd(1:tK), clNoNInd(1:tK)];
    clNoNNames  = pS(i).G.Nodes.Name(clTop10Ind(:,2));
    clIndvNames = pS(i).G.Nodes.Name(clTop10Ind(:,1));
    dgTop10Ind = [dgIndvInd(1:tK), dgNoNInd(1:tK)];
    dgNoNNames = pS(i).G.Nodes.Name(dgTop10Ind(:,2));
    dgIndvNames = pS(i).G.Nodes.Name(dgTop10Ind(:,1));
    efTop10Ind = [efIndvInd(1:tK), efNoNInd(1:tK)];
    efNoNNames = pS(i).G.Nodes.Name(efTop10Ind(:,2));
    efIndvNames = pS(i).G.Nodes.Name(efTop10Ind(:,1));


    % Alluvial plots (Figure 4) ===========================================
    scaleTF = string(pS(i).G.Nodes.scale{1}) == ["Urban", "Commuter", "National"];
    scaleColi = scaleCols(scaleTF,:);
    alluvOpts.patchColors = patchColors{scaleTF};

    if ~tKallTF
        fAbc = figure;
        [~, drBc(i), dI2Cbc, dC2Ibc] = alluvTopk(pS(i).G.Nodes.Name(bcIndvInd), pS(i).G.Nodes.Name(bcNoNInd), tK, alluvOpts);
        title(sysNames(i) + ", Betweenness", Interpreter="latex")
        fArwbc = figure;
        [~, drRwbc(i), dI2Crwbc, dC2Irwbc] = alluvTopk(pS(i).G.Nodes.Name(rwbcIndvInd), pS(i).G.Nodes.Name(rwbcNoNInd), tK, alluvOpts);
        title(sysNames(i) + ", Passenger-Weighted Betweenness", Interpreter="latex")
        fAdc = figure;
        [~, drDc(i), dI2Cdc, dC2Idc] = alluvTopk(pS(i).G.Nodes.Name(dgIndvInd), pS(i).G.Nodes.Name(dgNoNInd), tK, alluvOpts);
        title(sysNames(i) + ", Degree", Interpreter="latex")
        fAcc = figure;
        [~, drCc(i), dI2Ccc, dC2Icc] = alluvTopk(pS(i).G.Nodes.Name(clIndvInd), pS(i).G.Nodes.Name(clNoNInd), tK, alluvOpts);
        title(sysNames(i) + ", Closeness", Interpreter="latex")
        fAed = figure;
        title("By $\delta E$ node deletion")
        [~, drEd(i), dI2Ced, dC2Ied] = alluvTopk(pS(i).G.Nodes.Name(efIndvInd), pS(i).G.Nodes.Name(efNoNInd), tK, alluvOpts);
        title(sysNames(i) + ", $\delta$ Efficiency", Interpreter="latex")
        fArwbc.Position = [ 772, -26, 767, 1100];
        if ~alluvTF
            close(fAbc, fArwbc, fAdc, fAcc, fAed)
            % else
            %     close(fAdc, fAcc, fAed, fAbc) % Uncomment to produce only RWBC alluvial plots
        end
    end

    dChg_arr_bc{i,1} = dC2Ibc; dChg_arr_bc{i,2} = dI2Cbc; % Store results
    dChg_arr_rwbc{i,1} = dC2Irwbc; dChg_arr_rwbc{i,2} = dI2Crwbc;
    dChg_arr_dc{i,1} = dC2Idc; dChg_arr_dc{i,2} = dI2Cdc;
    dChg_arr_cc{i,1} = dC2Icc; dChg_arr_cc{i,2} = dI2Ccc;
    dChg_arr_ed{i,1} = dC2Ied; dChg_arr_ed{i,2} = dI2Ced;


    % Rank-shift analysis =================================================
    if ~tKallTF
        dI2Cbc = dI2Cbc(1:tK); dC2Ibc = dC2Ibc(1:tK); % Trim to top K ranks
        dI2Crwbc = dI2Crwbc(1:tK); dC2Irwbc = dC2Irwbc(1:tK);
        dI2Ccc = dI2Ccc(1:tK); dC2Icc = dC2Icc(1:tK);
        dI2Cdc = dI2Cdc(1:tK); dC2Idc = dC2Idc(1:tK);
        dI2Ced = dI2Ced(1:tK); dC2Ied = dC2Ied(1:tK);
    end

    [minDiffbc, maxDiffbc] = deal(dI2Cbc + dC2Ibc); % Preallocate
    [minDiffrwbc, maxDiffrwbc] = deal(dI2Crwbc + dC2Irwbc);
    [minDiffcc, maxDiffcc] = deal(dI2Ccc + dC2Icc);
    [minDiffdc, maxDiffdc] = deal(dI2Cdc + dC2Idc);
    [minDiffed, maxDiffed] = deal(dI2Ced + dC2Ied);

    minDiffbc(dI2Cbc.*dC2Ibc < 0) = min(dI2Cbc(dI2Cbc.*dC2Ibc < 0), dC2Ibc(dI2Cbc.*dC2Ibc < 0));
    maxDiffbc(dI2Cbc.*dC2Ibc < 0) = max(dI2Cbc(dI2Cbc.*dC2Ibc < 0), dC2Ibc(dI2Cbc.*dC2Ibc < 0));
    minDiffrwbc(dI2Crwbc.*dC2Irwbc < 0) = min(dI2Crwbc(dI2Crwbc.*dC2Irwbc < 0), dC2Irwbc(dI2Crwbc.*dC2Irwbc < 0));
    maxDiffrwbc(dI2Crwbc.*dC2Irwbc < 0) = max(dI2Crwbc(dI2Crwbc.*dC2Irwbc < 0), dC2Irwbc(dI2Crwbc.*dC2Irwbc < 0));
    minDiffcc(dI2Ccc.*dC2Icc < 0) = min(dI2Ccc(dI2Ccc.*dC2Icc < 0), dC2Icc(dI2Ccc.*dC2Icc < 0));
    maxDiffcc(dI2Ccc.*dC2Icc < 0) = max(dI2Ccc(dI2Ccc.*dC2Icc < 0), dC2Icc(dI2Ccc.*dC2Icc < 0));
    minDiffdc(dI2Cdc.*dC2Idc < 0) = min(dI2Cdc(dI2Cdc.*dC2Idc < 0), dC2Idc(dI2Cdc.*dC2Idc < 0));
    maxDiffdc(dI2Cdc.*dC2Idc < 0) = max(dI2Cdc(dI2Cdc.*dC2Idc < 0), dC2Idc(dI2Cdc.*dC2Idc < 0));
    minDiffed(dI2Ced.*dC2Ied < 0) = min(dI2Ced(dI2Ced.*dC2Ied < 0), dC2Ied(dI2Ced.*dC2Ied < 0));
    maxDiffed(dI2Ced.*dC2Ied < 0) = max(dI2Ced(dI2Ced.*dC2Ied < 0), dC2Ied(dI2Ced.*dC2Ied < 0));

    mindbc = min(minDiffbc); maxdbc = max(maxDiffbc);
    mindrwbc = min(minDiffrwbc); maxdrwbc = max(maxDiffrwbc);
    minddc = min(minDiffdc); maxddc = max(maxDiffdc);
    mindcc = min(minDiffcc); maxdcc = max(maxDiffcc);
    minded = min(minDiffed); maxded = max(maxDiffed);


    % Rank-shift plots (Figure 3) =========================================
    ytickLbls = "#" + string([1, floor(tK/2), tK]');
    minLimrwbc = min([dI2Crwbc; dC2Irwbc]);
    maxLimrwbc = max([dI2Crwbc; dC2Irwbc]);
    xtrwbc = [minLimrwbc, maxLimrwbc];

    if ~tKallTF
        ti = nexttile(tl1, tileOrder(i));
        bhi = barh(1:tK, [dI2Crwbc, dC2Irwbc]); hold on
        if prod(xtrwbc) > 0
            xticks(sort([0, xtrwbc(abs(xtrwbc) == max(abs(xtrwbc)))]))
        else
            xticks(sort([xtrwbc, 0]))
        end
        yticks(sort([1, floor(tK/2), tK]))
        yticklabels(ytickLbls)
        set(gca, 'YDir','reverse')
        set(bhi, "FaceAlpha", 0.7, "EdgeAlpha", 0)
        bhi(1).FaceColor = scaleColi;
        bhi(2).FaceColor = Ccol;
        ti.YAxisLocation = "right";

        title(ti, sysNames(i), Interpreter="latex")
        if i == numel(S)
            legend(["$\delta R_{I \rightarrow C} = R^{\sigma}_{I} \left ( V^{\sigma,k}_{I} \right ) - R^{\sigma}_{C} \left ( V^{\sigma,k}_{I} \right )$",...
                "$\delta R_{C \rightarrow I} = R^{\sigma}_{C} \left ( V^{\sigma,k}_{C} \right ) - R^{\sigma}_{I} \left ( V^{\sigma,k}_{C} \right )$"],...
                Interpreter="latex", FontSize=14)
        end
    end

    % Simulate sequential failure and recovery ============================
    [effGsNoNbc, effGsIndvbc, effSsNoNbc, effSsIndvbc] = seqFailRecovNEC(G, S{i}, bcIndvNames, bcNoNNames, resOpts);
    [effGsNoNrwbc, effGsIndvrwbc, effSsNoNrwbc, effSsIndvrwbc] = seqFailRecovNEC(G, S{i}, rwbcIndvNames, rwbcNoNNames, resOpts);
    [effGsNoNdc, effGsIndvdc, effSsNoNdc, effSsIndvdc] = seqFailRecovNEC(G, S{i}, dgIndvNames, dgNoNNames, resOpts);
    [effGsNoNcc, effGsIndvcc, effSsNoNcc, effSsIndvcc] = seqFailRecovNEC(G, S{i}, clIndvNames, clNoNNames, resOpts);
    [effGsNoNed, effGsIndved, effSsNoNed, effSsIndved] = seqFailRecovNEC(G, S{i}, efIndvNames, efNoNNames, resOpts);

    rsEffGsbc = [effGsNoNbc; effGsIndvbc]; rsEffSsbc = [effSsNoNbc; effSsIndvbc]; % store results
    rsEffGsrwbc = [effGsNoNrwbc; effGsIndvrwbc]; rsEffSsrwbc = [effSsNoNrwbc; effSsIndvrwbc];
    rsEffGsdc = [effGsNoNdc; effGsIndvdc]; rsEffSsdc = [effSsNoNdc; effSsIndvdc];
    rsEffGscc = [effGsNoNcc; effGsIndvcc]; rsEffSscc = [effSsNoNcc; effSsIndvcc];
    rsEffGsed = [effGsNoNed; effGsIndved]; rsEffSsed = [effSsNoNed; effSsIndved];
    
    % Failure and recovery plots (Figure 5) ===============================
    if tKallTF
        GfrPlt = rsEffGsrwbc; SfrPlt = rsEffSsrwbc;
        if any(sysNames(i) == ["NYC Subway", "Metro-North", "NJ Transit", "SEPTA Regional"])
            LblVAi = "top"; LblHAi = "left";
        else
            LblVAi = "middle"; LblHAi = "left";
        end

        tj = nexttile(tl2, tileOrder(i));
        plot(GfrPlt(1,:), '-', LineWidth=2, Color=Ccol); hold on
        plot(GfrPlt(2,:), ':', LineWidth=2, Color=Ccol);
        plot(SfrPlt(1,:), "-", LineWidth=2, Color=scaleColi);
        plot(SfrPlt(2,:), ":", LineWidth=2, Color=scaleColi);
        xline((numel(GfrPlt(1,:))+1)/2, "-k", "$t_{f \textrm{-} r}$", Interpreter="latex", FontSize=textSz, LabelVerticalAlignment=LblVAi, LabelHorizontalAlignment=LblHAi)
        title(tj, sysNames(i), Interpreter="latex", FontSize=textSz)
        ylim([-0.02, 1.1])
        xlim([0, numel(effGsNoNbc)+1])
        tj.Box = "off";
        xline(tj,tj.XLim(2))
        yline(tj,tj.YLim(2))
        set(tj,'TickLength',[0.01 0])
        tj.XTick = [0, tj.XTick(end)/2, tj.XTick(end)];
        tj.YTick = [0 1];

        figure;
        plt1 = plot(GfrPlt(1,:), '-', LineWidth=3, Color=Ccol); hold on
        plt2 = plot(GfrPlt(2,:), ':', LineWidth=3, Color=Ccol);
        plt3 = plot(SfrPlt(1,:), "-", LineWidth=3, Color=scaleColi);
        plt4 = plot(SfrPlt(2,:), ":", LineWidth=3, Color=scaleColi);
        xline((numel(GfrPlt(1,:))+1)/2, "-k", "$t_{f \textrm{-} r}$", Interpreter="latex", FontSize=textSz, LabelVerticalAlignment=LblVAi, LabelHorizontalAlignment=LblHAi)
        xline(1, "-k", "$t_{0}$", Interpreter="latex", FontSize=textSz, LabelVerticalAlignment=LblVAi, LabelHorizontalAlignment=LblHAi)
        xline(numel(GfrPlt(1,:)), "-k", "$t_{r}$", Interpreter="latex", FontSize=textSz, LabelVerticalAlignment=LblVAi, LabelHorizontalAlignment="right")
        title(sysNames(i), Interpreter="latex", FontSize=textSz)
        ax = gca; ax.Box = "off";
        set(ax,'TickLength',[0.01 0]); ax.XTick = [0, ax.XTick(end)/2, ax.XTick(end)]; ax.YTick = [0, 0.5, 1];
        xlabel("Step ($t$)", Interpreter="latex", FontSize=textSz); ylabel("$\rho_{\textrm{inst}}$", Interpreter="latex", FontSize=textSz+4)
        if i == numel(S)
            legend([plt1, plt2, plt3, plt4], ["$G^{\sigma}_{C}(t)$", "$G^{\sigma}_{I}(t)$", "$S^{\sigma}_{C}(t)$", "$S^{\sigma}_{I}(t)$"], Interpreter="latex", FontSize=textSz)
        end
    end

    % Store results =======================================================
    E_arr_bc{i,1} = rsEffGsbc; E_arr_bc{i,2} = rsEffSsbc;
    E_arr_rwbc{i,1} = rsEffGsrwbc; E_arr_rwbc{i,2} = rsEffSsrwbc;
    E_arr_dc{i,1} = rsEffGsdc; E_arr_dc{i,2} = rsEffSsdc;
    E_arr_cc{i,1} = rsEffGscc; E_arr_cc{i,2} = rsEffSscc;
    E_arr_ed{i,1} = rsEffGsed; E_arr_ed{i,2} = rsEffSsed;
    
end

E_arr_all = {E_arr_dc, E_arr_cc, E_arr_ed, E_arr_bc, E_arr_rwbc};

%% Resilience quantification ==============================================

FRmatAll = calcMultiResFR(E_arr_all); % dim 3: DC, CC, ED, BC, RWBC

rRfun = @(xi,yi) cellfun(@(xj,yj) (sum(xj(1, ((numel(xj(1,:))+1)/2):end)) / sum(xj(2, ((numel(xj(2,:))+1)/2):end)))...
    /(sum(yj(2, ((numel(yj(2,:))+1)/2):end)) / sum(yj(1, ((numel(yj(1,:))+1)/2):end))), xi, yi); % lambda_r
rRbc = rRfun(E_arr_bc(:,1), E_arr_bc(:,2));
rRrwbc = rRfun(E_arr_rwbc(:,1), E_arr_rwbc(:,2));
rRdc = rRfun(E_arr_dc(:,1), E_arr_dc(:,2));
rRcc = rRfun(E_arr_cc(:,1), E_arr_cc(:,2));
rRed = rRfun(E_arr_ed(:,1), E_arr_ed(:,2));

eOrd = [9 1:8 10:15];
sysCats = categorical("$\quad$" + sysNames(eOrd));
sysCats = reordercats(sysCats, "$\quad$" + sysNames(eOrd));
fullViewTF = true;

if miscPlotsTF
    figure;
    rRbar = barh(sysCats, [rRrwbc(eOrd), rRbc(eOrd), rRdc(eOrd), rRcc(eOrd), rRed(eOrd)], Interpreter="latex", FaceAlpha=0.8); hold on
    set(rRbar, "FaceColor", "flat")
    for ii = 1:5; rRbar(ii).CData(10:15,:) = urbColsQuad(ii,:).*ones(6,3); end
    for ii = 1:5; rRbar(ii).CData(2:9,:) = regColsQuad(ii,:).*ones(8,3); end
    for ii = 1:5; rRbar(ii).CData(1,:) = intColsQuad(ii,:); end
    set(rRbar, "FaceAlpha", 0.7, "EdgeAlpha", 0.2)
    xline(1, "-", "$\bf{\beta}_{r}$", Interpreter="latex", Alpha=0.3, LabelVerticalAlignment="top", FontSize=14)
    rRbar(1).Parent.YTick = [];
    rRbar(5).FontSize = 12; box off
    ax = rRbar.Parent; xt = ax.XTick;
    xt(2:2:numel(xt)) = []; ax.XTick = xt;
    xlabel("$\lambda_{r} (G,S)$", Interpreter="latex", FontSize=20)
end

%% Analysis of results ====================================================

% FR: double of size [(num systems) X (G_global, G_local, S_global, S_local) X (num strategies/metrics) X (recov phase, fail phase)]
FRmat_rGpGrec = squeeze(FRmatAll(:,1,:,1)); % R-G + S-G + recovery
FRmat_rGpLrec = squeeze(FRmatAll(:,2,:,1)); % R-G + S-L + recovery
FRmat_rGpGfail = squeeze(FRmatAll(:,1,:,2)); % R-G + S-G + failure
FRmat_rGpLfail = squeeze(FRmatAll(:,2,:,2)); % R-G + S-L + failure

FRmat_rLpGrec = squeeze(FRmatAll(:,3,:,1)); % R-L + S-G + recovery
FRmat_rLpLrec = squeeze(FRmatAll(:,4,:,1)); % R-L + S-L + recovery
FRmat_rLpGfail = squeeze(FRmatAll(:,3,:,2)); % R-L + S-G + failure
FRmat_rLpLfail = squeeze(FRmatAll(:,4,:,2)); % R-L + S-L + failure

% Metrics that maximize R-G
[maxv_rGpGrec, maxi_rGpGrec] = max(FRmat_rGpGrec, [], 2); 
[maxv_rGpLrec, maxi_rGpLrec] = max(FRmat_rGpLrec, [], 2);
[maxv_rGpGfail, maxi_rGpGfail] = max(FRmat_rGpGfail, [], 2);
[maxv_rGpLfail, maxi_rGpLfail] = max(FRmat_rGpLfail, [], 2);

% Metrics that maximize R-L
[maxv_rLpGrec, maxi_rLpGrec] = max(FRmat_rLpGrec, [], 2);
[maxv_rLpLrec, maxi_rLpLrec] = max(FRmat_rLpLrec, [], 2);
[maxv_rLpGfail, maxi_rLpGfail] = max(FRmat_rLpGfail, [], 2);
[maxv_rLpLfail, maxi_rLpLfail] = max(FRmat_rLpLfail, [], 2);

% Compare efficacy of S-G vs. S-L for all metrics
pctIncrGrec_all = 100.*(FRmat_rGpGrec-FRmat_rGpLrec)./FRmat_rGpLrec;
pctIncrGfail_all = 100.*(FRmat_rGpGfail-FRmat_rGpLfail)./FRmat_rGpLfail;
pctDecrLrec_all = 100.*(FRmat_rLpGrec-FRmat_rLpLrec)./FRmat_rLpLrec;
pctDecrLfail_all = 100.*(FRmat_rLpGfail-FRmat_rLpLfail)./FRmat_rLpLfail;

% Recov: Percent increase in R-G: best S-G metric vs. corresponding S-L metric
% Intercity: 12.0421% | urban: 14.9194% | regional: 22.0422% 
pctIncrGrec = 100*(maxv_rGpGrec-FRmat_rGpLrec(sub2ind(size(FRmat_rGpLrec), ...
    1:size(FRmat_rGpLrec,1), maxi_rGpGrec'))')./FRmat_rGpLrec(sub2ind(size(FRmat_rGpLrec), 1:size(FRmat_rGpLrec,1), maxi_rGpGrec'))';
pctIncrGrec_scaleAvg = [mean(pctIncrGrec(10:15)), mean(pctIncrGrec(1:8)), pctIncrGrec(9)];

% Failure: Percent increase in R-G: best S-G metric vs. corresponding S-L metric
% Intercity: 23.1907% | urban: 141.0308% | regional: 101.7578% 
pctIncrGfail = 100*(maxv_rGpGfail-FRmat_rGpLfail(sub2ind(size(FRmat_rGpLfail), ...
    1:size(FRmat_rGpLfail,1), maxi_rGpGfail'))')./FRmat_rGpLfail(sub2ind(size(FRmat_rGpLfail), 1:size(FRmat_rGpLfail,1), maxi_rGpGfail'))';
pctIncrGfail_scaleAvg = [mean(pctIncrGfail(10:15)), mean(pctIncrGfail(1:8)), pctIncrGfail(9)];

% Benefit-cost ratios (Increase in R-G / decrease in S-L), i.e., BCR = dRG/dSL
incGr_decLr_rec = pctIncrGrec_all./pctDecrLrec_all;
incGr_decLr_fail = pctIncrGfail_all./pctDecrLfail_all;

avgByMetricRec = mean(abs(incGr_decLr_rec),1); % average R for each metric, recovery
avgByMetricFail = mean(abs(incGr_decLr_fail),1); % average R for each metric, failure

[maxPctIncGrec, miPctIncGrec] = max(pctIncrGrec_all,[],2); % Maximum increase in R (recovery)
abs(maxPctIncGrec./arrayfun(@(xi,yi) pctDecrLrec_all(xi,yi), (1:size(miPctIncGrec,1))', miPctIncGrec))

% Best performance ratios for scale-averaged R
% Values correspond to Figure 6 (bottom right panel)
% FAILURE: Urban, degree; Regional, degree; Intercity, degree
topTradeoffsFail = [(0.300065-0.212562)/(.81204-.837676), (0.263821-0.186738)/(0.805687-0.823586), (0.630281-0.511631)/(0.805689-0.84916)];
toprGincFail = 100.*[(0.300065-0.212562)/0.212562, (0.263821-0.186738)/0.186738, (0.630281-0.511631)/0.511631];

% RECOVERY: Urban, closeness; Regional, closeness; Intercity, betweenness
% Values correspond to Figure 6 (top right panel)
topTradeoffsRecov = [(0.92057-0.807415)/(0.392184-0.437369), (0.921091-0.815495)/(0.40374-0.437302), (0.842685-0.752115)/(0.419188-0.442668)];
toprGincRecov = 100.*[(0.92057-0.807415)/0.807415, (0.921091-0.815495)/0.815495, (0.842685-0.752115)/0.752115];

% closeness, recovery, intercity
(0.765084-0.759401)/(0.465006-0.469378)
% Closeness abs percent over avg, recov, urban
100*mean(FRmat_rGpGrec(10:15, 2))/mean(FRmat_rGpGrec(10:15, :), "all")

% BC and RWBC, recov
BCRfun = @(y2,y1,x2,x1) (y2-y1)/(x2-x1);
% Unweighted betweenness BCR
BCRfun(mean(FRmat_rGpGrec(10:15, 4)), mean(FRmat_rGpLrec(10:15, 4)), mean(FRmat_rLpGrec(10:15, 4)), mean(FRmat_rLpLrec(10:15, 4)))
% RWBC BCR
BCRfun(mean(FRmat_rGpGrec(10:15, 5)), mean(FRmat_rGpLrec(10:15, 5)), mean(FRmat_rLpGrec(10:15, 5)), mean(FRmat_rLpLrec(10:15, 5)))


%% More analysis of results ===============================================

nNodesS = cellfun(@(x) height(x.Nodes), S);
E_lbls_all = ["Degree", "Closeness", "Efficiency", "Betweenness", "Betweenness (weighted)"];

% Recovery phase
AaRatio = @(x) cellfun(@(y) sum(y(1, ((numel(y(1,:))+1)/2):end))/((numel(y(1,:))+1)/2), x);
AbRatio = @(x) cellfun(@(y) sum(y(2, ((numel(y(2,:))+1)/2):end))/((numel(y(1,:))+1)/2), x);
AaAbRatio = @(x) cellfun(@(y) (sum(y(1, ((numel(y(1,:))+1)/2):end)) / sum(y(2, ((numel(y(2,:))+1)/2):end))), x);

% R-G ratio for S-G/S-L
GsRatiorwbc = AaAbRatio(E_arr_rwbc(:,1)); SsRatiorwbc = AaAbRatio(E_arr_rwbc(:,2));
GsRatiobc = AaAbRatio(E_arr_bc(:,1)); SsRatiobc = AaAbRatio(E_arr_bc(:,2));
GsRatiocc = AaAbRatio(E_arr_cc(:,1)); SsRatiocc = AaAbRatio(E_arr_cc(:,2));
GsRatiodc = AaAbRatio(E_arr_dc(:,1)); SsRatiodc = AaAbRatio(E_arr_dc(:,2));
GsRatioed = AaAbRatio(E_arr_ed(:,1)); SsRatioed = AaAbRatio(E_arr_ed(:,2));

% Summary stats (recovery)
% AaAb: ratio of (R-G+S-G)/(R-G+S-L), (R-L+S-L)/(R-L+S-G)
AaAb_rwbc = [SsRatiorwbc(eOrd), GsRatiorwbc(eOrd)];
AaAb_bc = [SsRatiobc(eOrd), GsRatiobc(eOrd)];
AaAb_cc = [SsRatiocc(eOrd), GsRatiocc(eOrd)];
AaAb_dc = [SsRatiodc(eOrd), GsRatiodc(eOrd)];
AaAb_ed = [SsRatioed(eOrd), GsRatioed(eOrd)];
AaAb_all = cat(3, AaAb_dc, AaAb_cc, AaAb_ed, AaAb_bc, AaAb_rwbc);
[~, miAaAb] = max(AaAb_all, [], 3);
GMaxStrats = [E_lbls_all(miAaAb(:,2)); string(sysCats)];

AaBb_max = cell2mat(arrayfun(@(ii) AaAb_all(ii,:,miAaAb(ii,2)), (1:size(AaAb_all,1))', 'UniformOutput', false));
AaAb_max_reg = mean(AaBb_max(2:9,:), 1);
AaAb_max_urb = mean(AaBb_max(10:15,:), 1);
AaAb_rwbc_reg = mean(AaAb_rwbc(2:9,:), 1);
AaAb_rwbc_urb = mean(AaAb_rwbc(10:15,:), 1);

% rho_r (rhor): integrated recov resilience normalized by timesteps
rhor_rwbc = [AaRatio(E_arr_rwbc(eOrd,2)), AaRatio(E_arr_rwbc(eOrd,1))];
rhor_bc = [AaRatio(E_arr_bc(eOrd,2)), AaRatio(E_arr_bc(eOrd,1))];
rhor_dc = [AaRatio(E_arr_dc(eOrd,2)), AaRatio(E_arr_dc(eOrd,1))];
rhor_cc = [AaRatio(E_arr_cc(eOrd,2)), AaRatio(E_arr_cc(eOrd,1))];
rhor_ed = [AaRatio(E_arr_ed(eOrd,2)), AaRatio(E_arr_ed(eOrd,1))];
rhor_all = cat(3, rhor_dc, rhor_cc, rhor_ed, rhor_bc, rhor_rwbc);
[~, miRhor] = max(rhor_all, [], 3);

rhorGs_max = cell2mat(arrayfun(@(ii) rhor_all(ii,:,miRhor(ii,2)), (1:size(rhor_all,1))', 'UniformOutput', false));
rhorSs_max = cell2mat(arrayfun(@(ii) rhor_all(ii,:,miRhor(ii,1)), (1:size(rhor_all,1))', 'UniformOutput', false));

%% Descriptive statistics & visualization =================================

% Degree distribution 
figure("Name", "Degree distributions") 
sp1 = subplot(4,4,1); hold on
histogram(degree(G), BinEdges=-0.5:10.5)
ylabel("Counts"); xlabel("Degree"); title("Northeast Corridor")
xlim([0,11])
for i = 1:numel(S)
    subplot(4,4,1+i); hold on
    histogram(degree(S{i}), BinEdges=-0.5:10.5)
    title(sysNames(i))
    xlim(sp1.XLim)
end

% Ridership distribution (Supplementary Figures)
figure("Name", "NEC ridership distribution") 
histogram(G.Nodes.rshp, 200)
ylabel("Counts"); xlabel("Annual ridership"); title("Northeast Corridor")

figure("Name", "Indv ridership distributions")
for i = 1:numel(S)
    subplot(4,4,i); hold on
    histogram(S{i}.Nodes.(rshpVarNames(i)), 200)
    title(sysNames(i))
end

% Summary stats table (Table 1)
ssT = table(Size=[15,9], VariableTypes={'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double'},...
    VariableNames=["System", "Scale", "Nodes", "Edges", "Loops", "Avg Degree", "Avg Shortest Path", "Clustering Coeff", "Modularity Index"]);
ssT.System = sysNames';
ssT.Scale = ["Regional", "Regional", "Regional", "Regional", "Regional", "Regional", "Regional", "Regional", ...
    "Intercity", "Urban", "Urban", "Urban", "Urban", "Urban", "Urban"]';
ssT.Nodes = cellfun(@(x) height(x.Nodes),S);
ssT.Edges = cellfun(@(x) height(x.Edges),S);
ssT.Loops = cellfun(@(x) numel(allcycles(x, MaxNumCycles=10000)), S);
ssT.("Avg Degree") = cellfun(@(x) mean(degree(x)), S);

DS = cellfun(@(x) distances(x), S, 'UniformOutput', false);
DmaskS = cellfun(@(x,y) triu(true(numnodes(y)), 1) & isfinite(x), DS, S, 'UniformOutput', false);
ssT.("Avg Shortest Path") = cellfun(@(x,y) mean(x(y)), DS, DmaskS);

%% Plotting: barH rank-shift averages (Figure 3.a)
% dChg_arr_rwbc: Col 1 = C2I, Col 2 = I2C

% Using RWBC as importance measure...
% Get top-10-ranked nodes for each RTN
dChgTop10 = cellfun(@(x,y) [y(1:10), x(1:10)], dChg_arr_dc(:,1), dChg_arr_dc(:,2),  'UniformOutput', false);
dChgTop10 = cat(3, dChgTop10{:});
% Normalize rank shift magnitudes by number of stations (nodes) in each RTN
dChgTop10norm = cellfun(@(x,y,z) [y(1:10), x(1:10)]./height(z.Nodes), dChg_arr_dc(:,1), dChg_arr_dc(:,2), S,  'UniformOutput', false);
dChgTop10norm = cat(3, dChgTop10norm{:});

dChgAvg = mean(dChgTop10, 3); % Get average rank shifts across all RTNs
dChgAvgUrb = mean(dChgTop10(:, :, 10:end), 3); % Get avg rank shifts for each of 3 scales (URTS/RCR/IPR)
dChgAvgReg = mean(dChgTop10(:, :, 1:8), 3);
dChgAvgNorm = mean(dChgTop10norm, 3);
dChgAvgUrbNorm = mean(dChgTop10norm(:, :, 10:end), 3);
dChgAvgRegNorm = mean(dChgTop10norm(:, :, 1:8), 3);

figure("Name", "Rank Change Summary")
bhi = barh(1:10, dChgAvgNorm); hold on
ax = gca; xticks(ax.XTick(2:2:end))
yticklabels("#" + string(ax.YTick)); set(gca, 'YDir','reverse'); ax.YAxisLocation = "left";
set(bhi, "FaceAlpha", 0.6, "EdgeAlpha", 0)
bhi(1).FaceColor = hex2rgb("#B8098F"); bhi(2).FaceColor = Ccol;
plot([dChgAvgUrbNorm(:,1), dChgAvgRegNorm(:,1)]', ((1:10)'-[0.15,0.15])', Color=[hex2rgb("#B8098F"),0.5])
plot([dChgAvgUrbNorm(:,2), dChgAvgRegNorm(:,2)]', ((1:10)'+[0.15,0.15])', Color=[Ccol,0.7])
scatter(dChgAvgUrbNorm, (1:10)'+[-0.15,0.15], 35, "Filled", MarkerEdgeColor=scaleCols(1,:), MarkerFaceColor=scaleCols(1,:), MarkerFaceAlpha=0.7, MarkerEdgeAlpha=0)
scatter(dChgAvgRegNorm, (1:10)'+[-0.15,0.15], 35, "Filled", MarkerEdgeColor=scaleCols(2,:), MarkerFaceColor=scaleCols(2,:), MarkerFaceAlpha=0.7, MarkerEdgeAlpha=0)
xlim([-0.62,0.13])
ax.XTick = [ax.XTick, 0.1];
ax.XTickLabel = {'-50%'; '-30%'; '-10%'; '+10%'};

dChgAvgAll = 100.*mean(dChgAvgNorm, 1); % Cited in body of paper
% Avg global-->local rank shift: -29.97%
% Avg local-->global rank shift: -8.00%


%% Displacement plots
scCols = [hex2rgb("#28AD29").*ones(8,3); hex2rgb("#4542B5"); hex2rgb("#AD2828").*ones(6,3)];

% One metric plot (abs and change, 2x3, rwbc + closeness, FAIL + RECOV ) 
% =========================================================================
fcbd2 = figure("Name", "Displacement scatter plot, cls/btw/diff, fail + recov");
pIdx = [2,5]; FRmatCcRwbc = zeros(15,4,2,2); % system X metric X scale X phase
FRDiffCcRwbc = zeros(15,2,2,2);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact')
axList = gobjects(2,3);
mStyles = ["+", "o", "square", "^", "diamond"];

for ii = 1:2
    GsNoNR = FRmatAll(:, 1, pIdx(ii), 1); GsNoNF = FRmatAll(:, 1, pIdx(ii), 2);
    GsIndvR = FRmatAll(:, 2, pIdx(ii), 1); GsIndvF = FRmatAll(:, 2, pIdx(ii), 2);
    SsNoNR = FRmatAll(:, 3, pIdx(ii), 1); SsNoNF = FRmatAll(:, 3, pIdx(ii), 2);
    SsIndvR = FRmatAll(:, 4, pIdx(ii), 1); SsIndvF = FRmatAll(:, 4, pIdx(ii), 2);
    FRmatCcRwbc = FRmatAll(:,:,pIdx,:);

    spii = nexttile(ii); hold on 
    axLbls = ["$\rho_{\textrm{r }}(S)$", ""];
    if ii == 1, axLbls(2) = "$\rho_{\textrm{r }}(G)$"; end
    displacementPlot(SsNoNR, GsNoNR, SsIndvR, GsIndvR, scCols, E_lbls_all(pIdx(ii)) + " , Recovery", axLbls, mStyles(pIdx(ii)))
    
    if ii > 1
        sc3 = scatter(1000, 1000, 40, [0.3 0.3 0.3], MarkerFaceAlpha=0.8, LineWidth=1); sc4 = scatter(1000, 1000, 50, [0.3 0.3 0.3], "filled", MarkerFaceAlpha=0.8);
        l1 = legend([spii.Children(5), spii.Children(17), spii.Children(11), sc4, sc3], [ "Urban", "Regional", "Intercity", "Global $(G)$", "Local ($S$)"], Interpreter="latex", Location="southeast", FontSize=12);
        l1.NumColumns = 2; l1.Title.String = "$\quad \quad$ System scale $\quad \quad$    Prioritization strategy $\quad$";
        spii.YTickLabel = [];
    end
    axis(spii,'tight'); squareAxes(spii); daspect([1 1 1])
    xlim([0.2, 0.5]); ylim([0.5, 1])
    axList(ii,1) = spii;

    spii = nexttile(ii+3); hold on 
    axLbls = ["$1 - \rho_{\textrm{f }}(S)$", ""];
    if ii == 1, axLbls(2) = "$1 - \rho_{\textrm{f }}(G)$"; end
    displacementPlot(SsNoNF, GsNoNF, SsIndvF, GsIndvF, scCols, E_lbls_all(pIdx(ii)) + " , Failure", axLbls, mStyles(pIdx(ii)))
    axis(spii,'tight'); squareAxes(spii); daspect([1 1 1])
    xlim([0.58, 1]); ylim([0, 0.7])
    if ii > 1, spii.YTickLabel = []; end
    axList(ii,2) = spii;
end

lStyles = ["-", "-", "-", "-", "-"];
scaleInds = {1:8, 10:15, 9};
scaleStyles = ["-", ":", "--"];
stratCols = cell2mat(arrayfun(@(x) hex2rgb(x), "#" + ["F54927"; "BA7623"; "2399BA"; "8323BA"; "92BA23"], UniformOutput=false));
pIdx = [2, 4, 5];  pltArr = gobjects(numel(pIdx),2,3); % CC, BC, RWBC
axArr = gobjects(2,1); axArr(1) = nexttile(3); axArr(2) = nexttile(6);

for ii = 1:numel(pIdx)
    pii = pIdx(ii);
    for ik = 1:2
        nexttile(3*ik); hold on
        if ii == 1 && ik == 2
            subAx6 = axes(Position=[.715 .055 .11 .11], YAxisLocation="right", XAxisLocation="top");
            hold on; box on; daspect([1 1 1])
        end
        for ij = 1:numel(scaleInds)
            if isscalar(scaleInds{ij}) && scaleInds{ij} == 9 && ik == 2
                tgtAx = subAx6;
            else
                tgtAx = axArr(ik);
            end
            pltArr(ii,ik,ij) = plot(tgtAx, mean(FRmatAll(scaleInds{ij}, 3:4, pii, ik), 1)', mean(FRmatAll(scaleInds{ij}, 1:2, pii, ik), 1)',...
                lStyles(pii), LineWidth=5.3, Color=[scCols(scaleInds{ij}(1),:) 0.2]);
            scatter(tgtAx, mean(FRmatAll(scaleInds{ij}, 3, pii, ik)), mean(FRmatAll(scaleInds{ij}, 1, pii, ik)), 70, scCols(scaleInds{ij}(1),:), mStyles(pii), "filled", MarkerFaceAlpha=0.8)
            scatter(tgtAx, mean(FRmatAll(scaleInds{ij}, 4, pii, ik)), mean(FRmatAll(scaleInds{ij}, 2, pii, ik)), 55, scCols(scaleInds{ij}(1),:), mStyles(pii), MarkerFaceAlpha=0.8, LineWidth=1)
        end
        axes(axArr(ik)); daspect([1 1 1]);  %#ok<LAXES>
    end
end

xlim(axArr(1), [0.2808 0.4831]); ylim(axArr(1), [0.7201 0.9651])
subAx6.YTick = [0.57, 0.61]; subAx6.XTick = [0.73, 0.77];
axLbls = ["$\hat{\rho}_{\textrm{r }}(S)$", "$\hat{\rho}_{\textrm{r }}(G)$";...
    "$1 - \hat{\rho}_{\textrm{f }}(S)$", "$1 - \hat{\rho}_{\textrm{f }}(G)$"];
xlabel(axArr(1), axLbls(1,1), Interpreter="latex", FontSize=14)
xlabel(axArr(2), axLbls(2,1), Interpreter="latex", FontSize=14)
ylabel(axArr(1), axLbls(1,2), Interpreter="latex", FontSize=14)
ylabel(axArr(2), axLbls(2,2), Interpreter="latex", FontSize=14)

axes(axArr(1)); scArr = gobjects(3,1);
axLims = [axArr(1).XLim; axArr(1).YLim];
for jj = 1:numel(pIdx)
    pjj = pIdx(jj);
    scArr(jj) = scatter(axArr(1), 1000, 1000, 70, [0.3 0.3 0.3], mStyles(pjj), "filled", MarkerFaceAlpha=0.8);
end
axArr(1).XLim = axLims(1,:); axArr(1).YLim = axLims(2,:);
l3 = legend(scArr, ["Closeness", "Betweenness", "$\omega$-Betweenness"], Interpreter="latex", Location="southwest", FontSize=12);
l3.Title.String = "Metric";

axList = reshape(axList,[6,1]); axList(5:6) = axArr; pause(3)
for ij = 1:6
    if string(class(axList(ij))) == "matlab.graphics.axis.Axes"
        axList(ij).XTick = axList(ij).XTick(1:2:numel(axList(ij).XTick)); 
        axList(ij).YTick = axList(ij).YTick(1:2:numel(axList(ij).YTick)); 
    end
end


% Version with urban and regional separated ===============================
if miscPlotsTF
    figure("Name", "Displacement scatter plots, separated"); tiledlayout(2,5)
    E_arr_all = {E_arr_dc, E_arr_cc, E_arr_ed, E_arr_bc, E_arr_rwbc};
    E_lbls_all = ["Degree", "Closeness", "Efficiency", "Betweenness", "Betweenness (P-weighted)"];
    for ii = 1:5
        Eaii = E_arr_all{ii};
        GsNoNR = cellfun(@(xj) sum(xj(1, ((numel(xj(1,:))+1)/2):end)), Eaii(:,1));
        GsIndvR = cellfun(@(xj) sum(xj(2, ((numel(xj(2,:))+1)/2):end)), Eaii(:,1));
        SsNoNR = cellfun(@(yj) sum(yj(1, ((numel(yj(1,:))+1)/2):end)), Eaii(:,2));
        SsIndvR = cellfun(@(yj) sum(yj(2, ((numel(yj(2,:))+1)/2):end)), Eaii(:,2));

        spii = nexttile(ii); hold on; %grid on
        for ki = 1:9
            plot([SsNoNR(ki)./nNodesS(ki), SsIndvR(ki)./nNodesS(ki)]', [GsNoNR(ki)./nNodesS(ki), GsIndvR(ki)./nNodesS(ki)]', LineWidth=5.3, Color=[scCols(ki,:) 0.2])
        end
        sc2i = scatter(SsIndvR(1:9)./nNodesS(1:9), GsIndvR(1:9)./nNodesS(1:9), 40, scCols(1:9,:), MarkerFaceAlpha=0.8, LineWidth=1);
        sc1i = scatter(SsNoNR(1:9)./nNodesS(1:9), GsNoNR(1:9)./nNodesS(1:9), 50, scCols(1:9,:), "filled", MarkerFaceAlpha=0.8);

        spij = nexttile(ii+5); hold on; %grid on
        for ki = 10:15
            plot([SsNoNR(ki)./nNodesS(ki), SsIndvR(ki)./nNodesS(ki)]', [GsNoNR(ki)./nNodesS(ki), GsIndvR(ki)./nNodesS(ki)]', LineWidth=5.3, Color=[scCols(ki,:) 0.2])
        end
        sc2j = scatter(SsIndvR(10:15)./nNodesS(10:15), GsIndvR(10:15)./nNodesS(10:15), 40, scCols(10:15,:), MarkerFaceAlpha=0.8, LineWidth=1);
        sc1j = scatter(SsNoNR(10:15)./nNodesS(10:15), GsNoNR(10:15)./nNodesS(10:15), 50, scCols(10:15,:), "filled", MarkerFaceAlpha=0.8);
        if ii == 1
            ylabel([spii, spij], "$\rho_{\textrm{r}}(G)$", Interpreter="latex", FontSize=14)
        else
            spii.YTickLabel = []; spij.YTickLabel = [];
            spii.XTickLabel = []; spij.XTickLabel = [];
        end
        if ii == 3; xlabel([spii, spij], "$\rho_{\textrm{r}}(S)$", Interpreter="latex", FontSize=14); end
        title(spii, E_lbls_all(ii), Interpreter="latex")
        xlim([spii, spij], [0.2, 0.5])
        ylim([spii, spij], [0.5, 1.1])
        daspect(spii, [1 1 1]); daspect(spij, [1 1 1]);
    end
end

%% Composite resilience plot, quad separated (Figure 5.a) ==================
% Plots average efficiency difference under strategies (S-G - S-L)
% Plots percentile confidence intervals for URTS and RCR (not for single IPR)

figure("Name", "NoN vs. Indv resilience difference, quad");
for kk = 1:4; subplot(2,2,kk); hold on
    ylk = yline(0, '-k', Alpha=0.3);
    xline(0.5, '-k',"$t_{f \textrm{-} r}$", Interpreter="latex", FontSize=14, Alpha=0.3)
    ylim([-0.5, 0.6])
end
sp1 = subplot(2,2,1); ylabel("$\rho(S_{C}) - \rho(S_{I})$", Interpreter="latex", FontSize=14); 
sp2 = subplot(2,2,2); 
sp3 = subplot(2,2,3); ylabel("$\rho(G_{C}) - \rho(G_{I})$", Interpreter="latex", FontSize=14); xlabel("Step", Interpreter="latex", FontSize=14)
sp4 = subplot(2,2,4); xlabel("Step", Interpreter="latex", FontSize=14);  

subplot(2,2,1); hold on
[~, ~, cmapr] = fanChart(xlq, (vertcat(EDqUrb{:,2}))', 'mean');
plot(xlq, EDqUrbSsAvg, "-", Color=hex2rgb("#7D0000"), LineWidth=2)
colormap(gca, flipud(cmapr)); cbr = colorbar; clim([5,95]); cbr.Ticks = (10:10:90)';
cbr.Location = "south"; cbr.Label.String = "Percentile"; cbr.Label.Interpreter = "latex";
cbr.TickLabelInterpreter = "tex"; cbr.Label.FontSize = 11;

subplot(2,2,3); hold on
fanChart(xlq, (vertcat(EDqUrb{:,1}))', 'mean')
plot(xlq, EDqUrbGsAvg, "-", Color=hex2rgb("#7D0000"), LineWidth=2)

subplot(2,2,2); hold on
[~, ~, cmapg] = fanChart(xlq, (vertcat(EDqCom{:,2}))', 'mean', [], 'alpha', 0.2, 'colormap', {'shadesOfColor', [0 .8 0]});
plot(xlq, EDq{9,2}, "-", Color=hex2rgb("#0991B8"), LineWidth=2)
plot(xlq, EDqComSsAvg, "-", Color=hex2rgb("#006316"), LineWidth=2)
colormap(gca, flipud(cmapg)); cbg = colorbar; clim([5,95]); cbg.Ticks = (10:10:90)';
cbg.Location = "south"; cbg.Label.String = "Percentile"; cbg.Label.Interpreter = "latex";
cbg.TickLabelInterpreter = "tex"; cbg.Label.FontSize = 11;

subplot(2,2,4); hold on
fanChart(xlq, (vertcat(EDqCom{:,1}))', 'mean', [], 'alpha', 0.2, 'colormap', {'shadesOfColor', [0 .8 0]})
plot(xlq, EDq{9,1}, "-", Color=hex2rgb("#0991B8"), LineWidth=2)
plot(xlq, EDqComGsAvg, "-", Color=hex2rgb("#006316"), LineWidth=2)

sp1.XTick = []; sp2.XTick = []; sp3.XTick = []; sp4.XTick = [];
sp1.YTick = [-0.3, 0, 0.3]; sp3.YTick = [-0.3, 0, 0.3]; sp4.YTick = []; sp2.YTick = [];
subplot(sp1); ph1 = scatter(1000, 1000, MarkerFaceColor=scaleCols(1,:), MarkerEdgeColor=scaleCols(1,:), MarkerFaceAlpha=0.4, MarkerEdgeAlpha=0.4);
ph2 = scatter(1000, 1000, MarkerFaceColor=scaleCols(2,:), MarkerEdgeColor=scaleCols(2,:), MarkerFaceAlpha=0.4, MarkerEdgeAlpha=0.4);
ph3 = scatter(1000, 1000, MarkerFaceColor=scaleCols(3,:), MarkerEdgeColor=scaleCols(3,:), MarkerFaceAlpha=0.4, MarkerEdgeAlpha=0.4);
legend([ph1, ph2, ph3], ["Urban", "Regional", "Intercity"], Interpreter="latex", FontSize=14)
colorbar

%% Composite resilience plots (indv vs non)
% Alternate version of Figure 5.a
% Everything on 2 axes, no CIs, shows all data with bolded scale averages

if miscPlotsTF
    fc = figure("Name", "NoN vs. Indv resilience difference");
    for kk = 1:2; subplot(1,2,kk); hold on
        ylk = yline(0, '-k', Alpha=0.3);
        xline(0.5, '-k',"$t_{f \textrm{-} r}$", Interpreter="latex", FontSize=14, Alpha=0.3)
        ylim([-0.5, 0.6])
    end
    sp1 = subplot(1,2,1); ylabel("$\rho[S_{C}] - \rho[S_{I}]$", Interpreter="latex", FontSize=14); xlabel("Step", Interpreter="latex", FontSize=14)
    sp2 = subplot(1,2,2); xlabel("Step", Interpreter="latex", FontSize=14);  
    maxLen = max(cellfun(@(x) size(x,2), E_arr_rwbc),[],"all");
    EDiff = cellfun(@(x) x(1,:) - x(2,:), E_arr_rwbc, 'UniformOutput', false);
    xlin = cellfun(@(x) linspace(0, 1, numel(x)), EDiff, 'UniformOutput', false);
    xlq = linspace(0,1,1000);
    EDq = cellfun(@(x,y) interp1(x, y, xlq), xlin, EDiff, "UniformOutput", false);
    EDqUrb = EDq(10:15,:);
    EDqCom = EDq(1:8,:);
    EDqComGsAvg = mean(vertcat(EDqCom{:,1}),1);
    EDqComSsAvg = mean(vertcat(EDqCom{:,2}),1);
    EDqUrbGsAvg = mean(vertcat(EDqUrb{:,1}),1);
    EDqUrbSsAvg = mean(vertcat(EDqUrb{:,2}),1);

    for kj = 1:size(E_arr_rwbc,1)
        scaleTF = string(pS(kj).G.Nodes.scale{1}) == ["Urban", "Commuter", "National"];
        scaleColk = scaleCols(scaleTF,:);
        EGsk = E_arr_rwbc{kj,1};
        ESsk = E_arr_rwbc{kj,2};
        ESskDiff = ESsk(1,:) - ESsk(2,:); % NoN/global minus indv/local
        EGskDiff = EGsk(1,:) - EGsk(2,:); % NoN/global minus indv/local
        nk = numel(ESskDiff);
        xk = linspace(0, 1, nk);
        xkq = linspace(0,1,1000);
        ESskDq = interp1(xk, ESskDiff, xkq);
        EGskDq = interp1(xk, EGskDiff, xkq);

        subplot(1,2,1); hold on
        pSsk = scatter(xk, ESskDiff, MarkerFaceColor=scaleColk, MarkerEdgeColor=scaleColk, MarkerFaceAlpha=0.05, MarkerEdgeAlpha=0.05);
        subplot(1,2,2); hold on
        pGsk = scatter(xk, EGskDiff, MarkerFaceColor=scaleColk, MarkerEdgeColor=scaleColk, MarkerFaceAlpha=0.05, MarkerEdgeAlpha=0.05);
    end

    subplot(1,2,1); plot(xlq, EDq{9,2}, "-", Color=pcIstart, LineWidth=2)
    plot(xlq, EDqComSsAvg, "-", Color=pcRstart, LineWidth=2)
    plot(xlq, EDqUrbSsAvg, "-", Color=pcUstart, LineWidth=2)
    subplot(1,2,2); plot(xlq, EDq{9,1}, "-", Color=pcIstart, LineWidth=2)
    plot(xlq, EDqComGsAvg, "-", Color=pcRstart, LineWidth=2)
    plot(xlq, EDqUrbGsAvg, "-", Color=pcUstart, LineWidth=2)

    sp1.XTick = []; sp2.XTick = [];
    sp1.YTick = [-0.3, 0, 0.3]; sp2.YTick = [-0.3, 0, 0.3];

    subplot(sp1); ph1 = scatter(1000, 1000, MarkerFaceColor=scaleCols(1,:), MarkerEdgeColor=scaleCols(1,:), MarkerFaceAlpha=0.4, MarkerEdgeAlpha=0.4);
    ph2 = scatter(1000, 1000, MarkerFaceColor=scaleCols(2,:), MarkerEdgeColor=scaleCols(2,:), MarkerFaceAlpha=0.4, MarkerEdgeAlpha=0.4);
    ph3 = scatter(1000, 1000, MarkerFaceColor=scaleCols(3,:), MarkerEdgeColor=scaleCols(3,:), MarkerFaceAlpha=0.4, MarkerEdgeAlpha=0.4);
    legend([ph1, ph2, ph3], ["Urban", "Regional", "Intercity"], Interpreter="latex", FontSize=14)
end



