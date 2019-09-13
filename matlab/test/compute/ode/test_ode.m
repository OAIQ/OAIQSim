t = [0,10];
lambda = -0.0001+10*1i;
ns = [10,50,250,1250,6250,31250];
nns = length(ns);

ye = exp(lambda*t(2));

err = zeros(nns,1);

for i=1:nns
    progressbar(i/nns);
    n = ns(i);
    y = test_ode_mex(real(lambda),imag(lambda),n,t);
    err(i) = norm(y-ye)/norm(ye);
end

plot(log(ns),log(err))