function u = eval(obj,X,XSize)
    % This function will evaluate the texture for the sample points
    % defined in the array X.
    % This should be updated once the GridData object is ready to
    % go!
    L = obj.L;   % Bounding box
    N = obj.N;
    if(numel(N)==1)
        N = N*ones(1,obj.dim);
    end
    if(nargin<2) % No array provided
        if(obj.dim==2)
            xtemp = linspace(L(1,1),L(1,2),N(1));
            ytemp = linspace(L(2,1),L(2,2),N(2));
            [xx,yy] = meshgrid(xtemp,ytemp);
            X = [xx(:),yy(:)];
            XSize = N;
        elseif(obj.dim==3)
            xtemp = linspace(L(1,1),L(1,2),N(1));
            ytemp = linspace(L(2,1),L(2,2),N(2));
            ztemp = linspace(L(3,1),L(3,2),N(3));
            [xx,yy,zz] = meshgrid(xtemp,ytemp,ztemp);
            X = [xx(:),yy(:),zz(:)];
            XSize = N;
        end
    end
    X = X';
    nEval = size(X,2);
    c_vec = obj.centers';
    if(numel(obj.b0)==1)
        b0_vec = obj.b0*ones(1,obj.K);
    elseif(numel(obj.b0)==obj.K)
        b0_vec = reshape(obj.b0,[1,obj.K]);
    else
        error('Incorrectly formatted b0! Must be 1-by-K or K-by-1');
    end
    if(numel(obj.rb2)==1)
        % Assume isotropic. I don't think I'm allowing this any more?
        rb2_vec = obj.rb2*ones(1,obj.dim*obj.K);
    elseif(numel(obj.rb2)==2)
        % Possibly anisotropic, but diagonal (2D)
        if(obj.gpu==1)
            rb2_vec = repmat([1/obj.rb2(1);0;1/obj.rb2(2)],[1,obj.K]);
            nParam  = 4;
        else
            rb2_vec = repmat(obj.rb2,[1,obj.K]);
        end
    elseif(numel(obj.rb2)==3)
        % Possibly anisotropic, but diagonal (3D)
        if(obj.gpu==1)
            rb2_vec = repmat([1/obj.rb2(1);0;0;1/obj.rb2(2);0;1/obj.rb2(3)],[1,obj.K]);
            nParam  = 7;
        else
            rb2_vec = repmat(obj.rb2,[1,obj.K]);
        end
    elseif(all(size(obj.rb2)==[obj.dim,obj.dim]))
        % Uniformly anisotropic (each lump has the same cov. mtx)
        if(obj.gpu==1)
            temp = inv(obj.rb2);
            mask = tril(true(size(obj.rb2)));
            temp = temp(mask); 
            rb2_vec = repmat(temp,[1,obj.K]);
            nParam = 1+obj.dim*(obj.dim+1)/2;
        else
            error('Anisotropic not available for CPU!');
        end
    elseif(all(size(obj.rb2)==[obj.dim,obj.dim*obj.K]))
        % Fully anisotropic (each lump has its own cov. mtx)
        if(obj.gpu==1)
            nParam = 1+obj.dim*(obj.dim+1)/2;
            rb2_vec = zeros(nParam-1,obj.K);
            mask = tril(true([obj.dim,obj.dim]));
            for i=1:obj.K
                temp = inv(obj.rb2(:,(obj.dim*(i-1)+1):(obj.dim*i)));
                temp = temp(mask);
                rb2_vec(:,i) = temp;
            end
        else
            error('Anisotropic not available for CPU!');
        end
    else
        disp(size(obj.rb2))
        error('Incorrectly formatted rb2!');
    end

    if(obj.gpu)
        theta_vec = [b0_vec;rb2_vec];
        u = lumpy_mex_gpu(uint32(obj.dim),uint32(obj.K),uint32(nEval),uint32(nParam),single(c_vec(:)),single(X(:)),single(obj.B0),single(theta_vec(:)));
    else
        %error('Need to update rb2_vec for CPU');
        u = lumpy_mex_cpu(uint32(obj.dim),uint32(obj.K),uint32(nEval),c_vec(:),X(:),obj.B0,b0_vec(:),rb2_vec(:));
    end
    if(nargin<2)  % No array provided, reshaping to locally generated grid size
        u = reshape(u,XSize);
    end
    if(nargin==3) % User supplied a grid size, reshape to it
        u = reshape(u,XSize);
    end
end


