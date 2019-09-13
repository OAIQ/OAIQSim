function result = ParseInputs(~,varargin)
    p = inputParser;
    p.CaseSensitive = true;
    default_support = RectSupport;
    
    ispos           = @(x) (x>0);
    isvaliddim      = @(x) (x==1||x==2||x==3);
    isvalidsupport  = @(x) isa(x,'RectSupport');
    
    addParameter(p,'support',default_support,isvalidsupport);
    addParameter(p,'grid_size',[],@isnumeric);
    addParameter(p,'L',[],@isnumeric);
    addParameter(p,'data_array',[],@isnumeric);
    
    p.parse(varargin{:});
    
    dim = p.Results.support.dim; % local value only
    result.support = p.Results.support;
    result.L   = p.Results.L;
    result.grid_size   = p.Results.grid_size;
    result.data_array   = p.Results.data_array;
  
    if(isempty(result.grid_size))
        result.grid_size = [64,64];
        if(isempty(result.data_array))     
            result.data_array = zeros(result.grid_size);
        end
    else
        if(length(result.grid_size)~=(dim+1))
            error('Incorrect size for grid_size.  Must be [Nx,Np], [Nx,Ny,Np] or [Nx,Ny,Nz,Np].')
        end
        if(isempty(result.data_array))
            result.data_array = zeros(result.grid_size);
        end
    end
    
    if(isempty(p.Results.L))
        result.L = repmat([0,1],[dim,1]);
    else
        if((size(result.L,1)~=dim)||size(result.L,2)~=2)
            error('Incorrect size for L');
        end
        if(any((result.L(:,2)>result.L(:,1))==0))
            error('Incorrect ordering for L; entries in right column should be larger than left column')
        end
    end
    
    if(isempty(result.data_array))
        result.data_array = zeros(result.grid_size);
    else
        sz = size(result.data_array);
        if(~(length(sz)==dim||length(sz)==(dim+1)))
            error('Incorrect size for data_array!')
        end
        
    end

       
   
end    
        