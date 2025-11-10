function FR = calcMultiResFR(sysPerfTS)

% OUT
% FR: double of size [(1,2,...,15 RTNs) X (R-G/S-G, R-G/S-L, R-L/S-G, R-L/S-L) X (DC, CC, DE, BC, RWBC) X (recovery, failure)]
%                    [    (system)      X   (global/local res & strat combo)   X   (importance metric)  X     (sim phase)    ]
P = sysPerfTS; % system performance cell array
numStrats = size(P,2);
numSystems = size(P{1},1);
FR = zeros(numSystems, 4, numStrats, 2);

for ii = 1:numStrats
    Pii = P{ii};
    GsNoNR = cellfun(@(xj) sum(xj(1, ((numel(xj(1,:))+1)/2):end))./((numel(xj(1,:))+1)/2), Pii(:,1));
    GsIndvR = cellfun(@(xj) sum(xj(2, ((numel(xj(2,:))+1)/2):end))./((numel(xj(2,:))+1)/2), Pii(:,1));
    SsNoNR = cellfun(@(yj) sum(yj(1, ((numel(yj(1,:))+1)/2):end))./((numel(yj(1,:))+1)/2), Pii(:,2));
    SsIndvR = cellfun(@(yj) sum(yj(2, ((numel(yj(2,:))+1)/2):end))./((numel(yj(2,:))+1)/2), Pii(:,2));

    GsNoNF = cellfun(@(xj) 1 - sum(xj(1, (1:(numel(xj(1,:))+1)/2 - 1)))./((numel(xj(1,:))+1)/2 - 1), Pii(:,1));
    GsIndvF = cellfun(@(xj) 1 - sum(xj(2, (1:(numel(xj(2,:))+1)/2 - 1)))./((numel(xj(2,:))+1)/2 - 1), Pii(:,1));
    SsNoNF = cellfun(@(yj) 1 - sum(yj(1, (1:(numel(yj(1,:))+1)/2 - 1)))./((numel(yj(1,:))+1)/2 - 1), Pii(:,2));
    SsIndvF = cellfun(@(yj) 1- sum(yj(2, (1:(numel(yj(2,:))+1)/2 - 1)))./((numel(yj(2,:))+1)/2 - 1), Pii(:,2));

    FR(:,:,ii,1) = [GsNoNR, GsIndvR, SsNoNR, SsIndvR];
    FR(:,:,ii,2) = [GsNoNF, GsIndvF, SsNoNF, SsIndvF];
end

end