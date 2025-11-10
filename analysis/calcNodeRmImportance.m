function P = calcNodeRmImportance(G, params, gFunc)

N = height(G.Nodes);
P = zeros(N,1);

for i = 1:N
    yprime = ones(N,1);
    yprime(i) = 0;
    P(i) = gFunc(yprime, params);
end

end