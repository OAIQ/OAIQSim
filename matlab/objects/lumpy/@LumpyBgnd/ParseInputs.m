function p = ParseInputs(~,varargin)
    p = inputParser;
    p.CaseSensitive = true;
    Kbar_default    = 100;
    b0_default      = 1;
    B0_default      = 0;
    cov_default     = 0.005;
    centers_default = [];
    if(~isempty(which('gpuDeviceCount')))
        gpu_default     = sign(gpuDeviceCount);
    else
        gpu_default = 0;
    end
    N_default       = 128;
    support_default = RectSupport;
    isnonneg        = @(x) (x>=0);
    ispos           = @(x) (x>0);
    isbinary        = @(x) (x==0||x==1);
    addParameter(p,'b0',b0_default,isnonneg);
    addParameter(p,'B0',B0_default,isnonneg);
    addParameter(p,'cov',cov_default,ispos);
    addParameter(p,'Kbar',Kbar_default,ispos);
    addParameter(p,'centers',centers_default,@isnumeric);
    addParameter(p,'gpu',gpu_default,isbinary);
    addParameter(p,'N',N_default,ispos);
    addParameter(p,'support',support_default,@(x) isa(x,'SupportSet'));
    addParameter(p,'pdf',false,@islogical);
    parse(p,varargin{:});
end