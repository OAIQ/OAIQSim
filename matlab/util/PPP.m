function Z = PPP(lambda,Nbar,dim,domain,ub,os)
    % - lambda is the intensity function - expected to be vectorized to
    % accept an Mxdim array and return an Mx1 vector.
    % - Nbar is the expected number of points 
    % - domain is a matrix of size 2xdim containing the bounds e.g. 
    % [0,1] -> samples the interval [0,1]
    % [0,1;0,1] -> samples the unit square 
    % [0,1;0,1;0,1] -> samples the unit cube, etc.  Only rectangular
    % domains considered for now.
    % - ub is an upper bound on lambda (needed for the rejection sampler)
    % - os is an 'oversample' factor to increase number of samples used by
    % rejection sampler, hence reducing number of lambda function evaluations.
    rng('shuffle');
    N = poissrnd(Nbar);
    enough_samples = 0;
    
    Z = [];
    while(enough_samples==0)
        Y = zeros(os*N,dim+1);
        for i=1:dim
            h = domain(i,2)-domain(i,1);
            Y(:,i) = domain(i,1)+h*rand(os*N,1);
        end
        Y(:,end) = rand(os*N,1);
        U = lambda(Y(:,1:dim));
        Z = [Z;Y(Y(:,dim+1)<U/ub,1:dim)];
        disp(['Number of samples is ',num2str(size(Z,1)),'/',num2str(N)]);
        if(size(Z,1)>=N)
            enough_samples = 1;
            Z = Z(1:N,:);
        end
    end
end