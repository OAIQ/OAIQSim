function p = parseinputs(~,varargin)
    p = inputParser;
    p.CaseSensitive = true;
    geom_default_2d    = [1,.1];
    geom_default_3d    = [1,1,0.1];
    pos_default_2d     = [1,0];
    pos_default_3d     = [1.5,0,0.5];

    pin_default     = [2e-1,1];  % Default is 
    blur_default    = 1e-2;   
    pinc_default    = 0;
    dim_default     = 2;
    coord_default_2d = 'cart';
    coord_default_3d = 'cyl';
    isvalidcoord = @(x) any(strcmp(x,{'cart','cyl','sph'}));

    addParameter(p,'dim',    dim_default,@isnumeric);
    p.KeepUnmatched = true;
    parse(p,varargin{:});
    
    if(p.Results.dim==2)
        addParameter(p,'detl',   geom_default_2d,@isnumeric);
        addParameter(p,'detpos',  pos_default_2d,@isnumeric);
        addParameter(p,'coordsys', coord_default_2d,isvalidcoord);
    elseif(p.Results.dim==3)
        addParameter(p,'detl',   geom_default_3d,@isnumeric);
        addParameter(p,'detpos',  pos_default_3d,@isnumeric);
        addParameter(p,'coordsys', coord_default_3d,isvalidcoord);
    end

    addParameter(p,'pinc',  pinc_default,@isnumeric);
    addParameter(p,'pinw',pin_default(1),@isnumeric);
    addParameter(p,'pind',pin_default(2),@isnumeric);
    addParameter(p,'blur',  blur_default,@isnumeric);
    parse(p,varargin{:}); 
end