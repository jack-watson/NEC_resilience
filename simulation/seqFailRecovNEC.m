function [effGsNoN, effGsIndv,effSsNoN, effSsIndv] = seqFailRecovNEC(G, S, indvNames, NoNNames, resOpts)

if nargin < 5
    resOpts = struct;
    resOpts.weightByRshpTF = false;
end

GsNoN = removeAndRestore(G, NoNNames);
GsIndv = removeAndRestore(G, indvNames);
SsNoN = removeAndRestore(S, NoNNames);
SsIndv = removeAndRestore(S, indvNames);

if resOpts.weightByRshpTF % If we're weighting efficiency via passenger flow then:
    GNames = string(G.Nodes.Name); % all node names in G and S
    SNames = string(S.Nodes.Name);

    nsteps = (numel(indvNames)*2)+1; % failure + recovery + 1
    ODindv = resOpts.ODMatrix.indv; % extract OD matrices from resOpts struct
    ODNoN = resOpts.ODMatrix.NoN;    

    [isOnSsIndv, isOnSsNoN] = deal(ones(size(ODindv,1), nsteps)); % preallocate
    [isOnGsIndv, isOnGsNoN] = deal(ones(size(ODNoN,1), nsteps));

    for i = 2:numel(indvNames)+1 % remove ones one at a time in sequence
        tormIi = string(indvNames(i-1)); % node names to remove in sequence
        tormNi = string(NoNNames(i-1));
        isOnSsIndv(:,i) = isOnSsIndv(:,i-1) - double(SNames == tormIi);
        isOnSsNoN(:,i)  = isOnSsNoN(:,i-1)  - double(SNames == tormNi);
        isOnGsNoN(:,i)  = isOnGsNoN(:,i-1)  - double(GNames == tormNi);
        isOnGsIndv(:,i) = isOnGsIndv(:,i-1) - double(GNames == tormIi);
    end

    numNodes = numel(indvNames); 
    for i = 1:numNodes
        toaddIi = string(indvNames(i)); % names of nodes to add back to networks
        toaddNi = string(NoNNames(i));
        isOnSsIndv(:,i+1+numNodes) = isOnSsIndv(:,i+numNodes) + double(SNames == toaddIi);
        isOnSsNoN(:,i+1+numNodes)  = isOnSsNoN(:,i+numNodes)  + double(SNames == toaddNi);
        isOnGsIndv(:,i+1+numNodes) = isOnGsIndv(:,i+numNodes) + double(GNames == toaddIi);
        isOnGsNoN(:,i+1+numNodes)  = isOnGsNoN(:,i+numNodes)  + double(GNames == toaddNi);
    end

end

paramsS = struct;
paramsG = struct;
paramsS.G = S;
paramsS.numNodes = height(paramsS.G.Nodes);
paramsS.isWeighted = true;
paramsG.G = G;
paramsG.numNodes = height(G.Nodes);
paramsG.isWeighted = true;
paramsS.resOpts = resOpts;
paramsG.resOpts = resOpts;
paramsG.resOpts.ODMatrix.indv = paramsG.resOpts.ODMatrix.NoN;
[effGsNoN, effGsIndv, effSsNoN, effSsIndv] = deal(zeros(size(GsNoN)));
[paramsGsN, paramsGsI] = deal(paramsG);
[paramsSsI, paramsSsN] = deal(paramsS);

if resOpts.weightByRshpTF
    paramsGsN.resOpts.ODMatrix.indv = ODNoN;
    paramsGsI.resOpts.ODMatrix.indv = ODNoN;
    paramsSsI.resOpts.ODMatrix.indv = ODindv;
    paramsSsN.resOpts.ODMatrix.indv = ODindv;
end

for j = 1:numel(effGsNoN)
    yG = ones(height(GsNoN{j}.Nodes),1);
    yS = ones(height(SsNoN{j}.Nodes),1);
    paramsGsN.G = GsNoN{j};
    paramsGsI.G = GsIndv{j};
    paramsSsN.G = SsNoN{j};
    paramsSsI.G = SsIndv{j};

    if resOpts.weightByRshpTF
        paramsGsN.resOpts.ODMatrix.activeNodeIndices = isOnGsNoN(:,j);
        paramsGsI.resOpts.ODMatrix.activeNodeIndices = isOnGsIndv(:,j);
        paramsSsN.resOpts.ODMatrix.activeNodeIndices = isOnSsNoN(:,j);
        paramsSsI.resOpts.ODMatrix.activeNodeIndices = isOnSsIndv(:,j);
    end

    effGsNoN(j) = networkEfficiency(yG, paramsGsN);
    effGsIndv(j) = networkEfficiency(yG, paramsGsI);
    effSsNoN(j) = networkEfficiency(yS, paramsSsN);
    effSsIndv(j) = networkEfficiency(yS, paramsSsI);
end

resNorm = @(x) x./(x(1));
effGsNoN = resNorm(effGsNoN);
effSsNoN = resNorm(effSsNoN);
effSsIndv = resNorm(effSsIndv);
effGsIndv = resNorm(effGsIndv);

end