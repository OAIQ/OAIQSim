


S = RectSupport([0,0.1;0,0.1;0,0.01]);

L = LumpyBgnd('S',S,'gpu',1);
L.cov = [1e-5;1e-5;1e-6];

K = 2.^(5:9);
N = 2.^(5:9);
nN = length(N);
nK = length(K);
times = zeros(nN,nK);

for i=1:nN
    for j=1:nK
        L.N = N(i);
        L.Kbar = K(j);
        tic
        u = L.eval;
        times(i,j) = toc;
        progressbar(i/nN,j/nK);
    end
end

