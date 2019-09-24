function randomize(obj,varargin)
    
    % Randomize the LumpyBgnd parameters according to the distributions
    % specified
    % The input which_params is a cell array of strings that specifies 


    %  Should start calling the object's support function unifRand function
    %  instead - the support set might not be rectangular 
    L = obj.L;  % Get bounding box;  L = [a,b;c,d] for [a,b]x[c,d] or L = [a,b;c,d;e,f] for 
    %           
    
    
    
    if(any(ismember(varargin,'shuffle')))
        disp('Shuffling RNG');
        rng('shuffle');
    end
    
    which_params = {}; 
    if(any(ismember(varargin,'all')))
        which_params = {'centers','K','b'};
    end
    if(any(ismember(varargin,'centers')))
        which_params{end+1} = 'centers';
    end
    if(any(ismember(varargin,'K')))
        which_params{end+1} = 'K';
    end
    if(any(ismember(varargin,'b')))
        which_params{end+1} = 'b';
    end
    
    if(nargin<2)
        which_params = {'centers','K','b'};
    end
    
    which_params
    
    if(any(ismember(which_params,'K')))
        k = obj.K_distr(obj.Kbar);
    else 
        k = obj.K;
    end
    
    
    if(any(ismember(which_params,'centers')))
        padx = obj.pad_factor*sqrt(max(obj.cov(1,:))/2);  %padding to improve stationarity
        pady = obj.pad_factor*sqrt(max(obj.cov(2,:))/2);

        minx = L(1,1)-padx;  
        maxx = L(1,2)+padx;
        miny = L(2,1)-pady;
        maxy = L(2,2)+pady;
        if(obj.dim==3)
            padz = obj.pad_factor*sqrt(max(obj.cov(3,:))/2);
            minz = L(3,1)-padz;
            maxz = L(3,2)+padz;
        end

        obj.centers = zeros(k,obj.dim);
        c_distr = obj.centers_distr;
        obj.centers(:,1) = minx + (maxx-minx).*c_distr(k);
        obj.centers(:,2) = miny + (maxy-miny).*c_distr(k);
        if(obj.dim==3)
            obj.centers(:,3) = minz + (maxz-minz).*c_distr(k);
        end
    end
    
    
    if(any(ismember(which_params,'b')))
        obj.b = obj.b_distr(k);
    end
    
end