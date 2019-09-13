function z = plus(x,y)
% Computes the sum of two lumpy backgrounds

    suppx = x.support;
    suppy = y.support;
    Kx    = x.K;
    Ky    = y.K;
    Covx  = x.cov;
    Covy  = y.cov;
    Nx    = x.N;
    Ny    = y.N; 

    if(~all(suppx.L(:)==suppy.L(:)))
        error('Incompatible supports!')
    end


    z = LumpyBgnd('support',suppx,'N',max(Nx,Ny));
    z.SetRandomLumps(0);

    z.centers = [x.centers;y.centers];
    
    if(size(Covx,2)==x.dim)
        % Uniform cov, need to repeat it
        Covx = repmat(Covx,[1,Kx]);
    end
    if(size(Covy,2)==y.dim)
        % Uniform cov, need to repeat it
        Covy = repmat(Covy,[1,Ky]);
    end
    
    % Need to use non-uniform covariance matrices since in theory the two
    % backgrounds can have different covs!
    z.cov = [Covx,Covy];  
    
    
    %z.K = Kx + Ky;

    if(numel(x.b0)==1)
        b0x = x.b0*ones(Kx,1);
    elseif(numel(x.b0)==Kx)
        b0x = x.b0(:);
    else
        error('Incorrectly formatted b0!');
    end

    if(numel(y.b0)==1)
        b0y = y.b0*ones(Ky,1);
    elseif(numel(y.b0)==Ky)
        b0y = y.b0(:);
    else
        error('Incorrectly formatted b0!');
    end


    z.b0 = [b0x;b0y];



end