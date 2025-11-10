function [fig, difRatio, dA2B, dB2A] = alluvTopk(lNames, rNames, k, opts)

A = lNames;
B = rNames;

% Find indices bi so that A(bi) = B
[~, bi] = ismember(B, A);

% Find indices ai so that B(ai) = A
[~, ai] = ismember(A, B);

uqlNames = unique([lNames(1:k); rNames(1:k)], "stable");
uqN = numel(uqlNames);
data = zeros(uqN);

AinBidx = bi(ismember(B, uqlNames));
BinAidx = ai(ismember(A, uqlNames));
lNamesk = regexprep(B(BinAidx), '_ \d+', '');
rNamesk = regexprep(A(AinBidx), '_ \d+', '');

lOverIdx = find(BinAidx > k);
rOverIdx = find(AinBidx > k);
adjPos = (1:numel(rOverIdx)) + k;

[~,lSortIdx] = sort(BinAidx(lOverIdx), "ascend");
[~,rSortIdx] = sort(AinBidx(rOverIdx), "ascend");

AinBpos = AinBidx;
BinApos = BinAidx;

BinApos(lOverIdx) = adjPos(lSortIdx);
AinBpos(rOverIdx(rSortIdx)) = adjPos;

for i = 1:numel(BinApos)
    data(AinBpos(i), i) = 1;
end

if nargin < 4 || ~isstruct(opts)
    opts = struct; 
end

fig = alluvialflow(data, lNamesk, rNamesk, sort(AinBidx), sort(BinAidx), opts);

dA2B = sort(AinBidx) - BinAidx; % change isolated to connected
dB2A = sort(BinAidx) - AinBidx; % change connected to isolated
difRatio = mean(dB2A(1:k))/mean(dA2B(1:k));

end