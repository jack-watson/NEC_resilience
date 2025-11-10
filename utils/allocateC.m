function [C_al, C_rem] = allocateC(B, C, Yt, Xt, params)

dmg = 1 - Yt; % damage magnitude is 1 minus node integrity
% U is maximum possible allocation for each damaged node s.t. for node i,
% yi(t+1) = 1.
U = max(0, ceil(dmg ./ params.alpha_c) - Xt); 

if sum(U) <= C % If hypothetical max alloc is within C budget
    C_rem = C - sum(U); % then remainder is C - sum(U)
    %C = sum(U);
    C_al = U;
    return
else
    C_rem = 0;
end

if all(U == 0)
    C_al = U;
    if size(C_al,1) < size(C_al, 2)
        C_al = C_al';
    end
    return
end

alpha = params.alpha_b;

% Define raw (fractional) allocation
% ...weird but typical PL behavior where diff condition returns TRUE but 
% mean condition returns FALSE depending on the order of magnitude of B
if all(diff(B) == 0) || all(B == mean(B)) 
    nToAlloc = sum(U > 0);
    allocIdx = find(U > 0);
    Craw = zeros(numel(U),1);
    Craw(allocIdx) = ones(nToAlloc,1).*(C/nToAlloc); % NEED TO ACCOUNT FOR U=0 CASES IN INITIAL RAW ALLOCATION TO GET ALLOCUNIT... LEFT OFF HERE B4 AUGUST
    allocUnit = ceil(max(Craw));
    C_rem = C;
    C_al = zeros(numel(Craw),1);
    %i = 0;
    while sum(C_al) < C && C_rem > 0 
        %i = i + 1;
        randInt = randi([1, nToAlloc]);
        randIdx = allocIdx(randInt);
        allocIdx(randInt) = [];
        nToAlloc = nToAlloc - 1;
        alloci = min([U(randIdx), allocUnit, C_rem]);
        C_al(randIdx) = C_al(randIdx) + alloci;
        C_rem = C_rem - alloci;
        %U = max(0, ceil(dmg ./ params.alpha_c) - Xt - C_al); 
        if nToAlloc == 0 && C_rem > 0 && any(C_al - U < 0)
            allocUnit = 1;
            nToAlloc = sum((C_al - U).*-1 > 0);
            allocIdx = find((C_al - U) < 0);
        end
    end
else
    if all(abs(B - mean(B)) < 0.00001)
        B = rescale(B, 0.001, 0.999);
    end

    Craw = C.*(B.^alpha)./sum(B.^alpha);
    
    % If we're allocating more resources to a node than it takes to fully
    % repair that node in a single timestep, i.e. if C(i)*alpha_c + Yt(i) >= 1 + alpha_c,
    % then reduce node i's allocation to its ceiling U(i) and distribute
    % the remaining excess C to the other nodes, excluding i and any other
    % previously over-allocated nodes from the calculation and using the
    % same allocation scheme (same alpha_b). This will likely result in
    % another over-allocation since we use the same alpha_b that resulted
    % in over-allocation in the initial iteration, so we continue this
    % process until either all C is allocated or all nodes reach their
    % allocation ceiling.
    UCDiff = Craw - U;
    toAlloc = true(numel(UCDiff),1);

    while any(UCDiff > 0) % if we've over-allocated any resources...
        ovalIdx = find(UCDiff > 0); % find over allocation indices
        toAlloc(ovalIdx) = false;
        Craw(ovalIdx) = U(ovalIdx);
        Covr = C - sum(Craw);
        Craw(toAlloc) = Covr.*(B(toAlloc).^alpha)./sum(B(toAlloc).^alpha);
        UCDiff = Craw - U;
    end

    C0 = floor(Craw);
    C0(C0 < 0) = 0;
    D = C - sum(C0);
    R = Craw - C0;
    C_al = C0;
    
    if ~all(R == 0)
        [~, iR] = maxk(R, D);
        allocUnit = ceil(D/numel(iR));
        for i = 1:numel(iR) % i = 1:D
            alloci = min([U(iR(i)), allocUnit, D]);
            C_al(iR(i)) = C_al(iR(i)) + alloci;
            D = D - alloci;
        end
    end

    UCDiff = C_al - U;
    if any(UCDiff > 0)
        ovalIdx = find(UCDiff > 0, 1);
        C_al(ovalIdx) = C_al(ovalIdx) - UCDiff(ovalIdx);
    end

    C_rem = C - sum(C_al);
end


end