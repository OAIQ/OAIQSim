function p = ParseInputs(~,varargin)
    p = inputParser;
    p.CaseSensitive = true;
    Kbar_default    = 100;
    %b_default       = 1;
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
    K_distr_default = @(Kb) poissrnd(Kb);
    centers_distr_default = @(K) rand(K,1);
    b_distr_default = @(K) ones(K,1);
    
    isnonneg        = @(x) (x>=0);
    ispos           = @(x) (x>0);
    isbinary        = @(x) (x==0||x==1);
    %addParameter(p,'b',b_default,isnonneg);
    addParameter(p,'B0',B0_default,isnonneg);
    addParameter(p,'cov',cov_default,ispos);
    addParameter(p,'Kbar',Kbar_default,ispos);
    addParameter(p,'centers',centers_default,@isnumeric);
    addParameter(p,'gpu',gpu_default,isbinary);
    addParameter(p,'N',N_default,ispos);
    addParameter(p,'support',support_default,@(x) isa(x,'SupportSet'));
    addParameter(p,'pdf',false,@islogical);
    addParameter(p,'K_distr',K_distr_default,@(x) isa(x,'function_handle'));
    addParameter(p,'centers_distr',centers_distr_default,@(x) isa(x,'function_handle'));
    addParameter(p,'b_distr',b_distr_default,@(x) isa(x,'function_handle'));
    parse(p,varargin{:});
end