% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:     RTELumpy.m
% Author:   Nick Henscheid
% Date:     6-2019
% Info:     Class definition file for the RTELumpy object, a ballistic RTE
%           solver for lumpy background objects.
% Inputs:   An object of type RTELumpy has the following optional inputs:
%           - 'Xi' a source function of LumpyBgnd type
%           - 'Mu' an attenuation function of LumpyBgnd type
%           - 'cm' a medium speed of light
%           - 'hmin' a minimum step length for the ballistic RTE solver 
%           - Defaults will be chosen for each if not specified by user.
% Notes:    - The dimension is fixed by the Xi and Mu functions (which must
%             obviously agree)
%           - The support sets for Xi and Mu must also agree.
% To Do:    - 
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------
classdef RTELumpy < handle
    properties (SetObservable = true)
       Xi LumpyBgnd;  % Lumpy object for source distribution
       Mu LumpyBgnd;  % Lumpy object for attenuation function
       cm;            % Medium light speed
       hmin;          % Minimum step length for integrator
    end
    
    properties (SetAccess = private)
       gpu;
    end
    
    properties (Dependent)
        dim;
    end
    
    methods
        function obj = RTELumpy(varargin)
            p         = obj.parseinputs(varargin{:});
            obj.Xi    = p.Results.Xi;
            obj.Mu    = p.Results.Mu;
            obj.cm    = p.Results.cm;
            obj.hmin  = p.Results.hmin;
            hasgpu = logical(gpuDeviceCount);
            if(hasgpu)
                obj.gpu   = p.Results.gpu;
            else
                disp('No GPU found.  Using CPU!');
                obj.gpu   = false;
            end
        end
        
        function setGPU(obj,gpu)
            hasgpu = logical(gpuDeviceCount);
            if(hasgpu)
                obj.gpu   = logical(gpu);
            else
                disp('No GPU found.  Using CPU!');
                obj.gpu   = false;
            end
        end
        
        function val = get.dim(obj)
            val = obj.Xi.dim;
        end
        
        function w = eval(obj,R,S,E)
           % Evaluates w for vectors R,S and E.  R,S and E must be the same
           % length, dim*neval-by-1.  length(w) = neval.
           % It is assumed that R and S are already unrolled, column-major,
           % from dim-by-neval (not sure why, should process this)
           if(isequal(size(R),size(S)))
               if(size(R,2)>1||size(S,2)>1||size(E,2)>1)
                   error('Wrong size inputs!')
               else 
                   Xi_in = obj.ExtractLumpyParam(obj.Xi);
                   Mu_in = obj.ExtractLumpyParam(obj.Mu);
                   neval = length(R)/obj.dim;
                   L_in  = obj.Xi.support.L'; L_in = L_in(:);% Row major
                   if(obj.gpu)
                       fprintf('Computing with GPU, neval = %i\n',neval);
                       w = rte_lumpy_mex_gpu(obj.dim,Xi_in,Mu_in,R,S,uint32(neval),obj.cm,obj.hmin,L_in);
                   else
                       fprintf('Computing with CPU, neval = %i\n',neval);
                       nstep = uint32(ceil(1/obj.hmin));
                       w = rte_lumpy_mex_cpu(obj.dim,Xi_in,Mu_in,R,S,uint32(neval),obj.cm,nstep,L_in);
                   end
               end
           else
               error('Inputs must be same size');
           end
        end
        
   
        function p = parseinputs(~,varargin)
            p = inputParser;
            p.CaseSensitive = true;
            Xi_default       = LumpyBgnd;
            Mu_default       = LumpyBgnd;
            cm_default      = 299792458;
            hmin_default   = 0.01;   % Step size for integrator
            gpu_default     = false;
            islumpy = @(l)isa(l,'LumpyBgnd');
            
            addParameter(p,'Xi',Xi_default,islumpy);
            addParameter(p,'Mu',Mu_default,islumpy);
            addParameter(p,'cm',cm_default,@isnumeric);
            addParameter(p,'hmin',hmin_default,@isnumeric);
            addParameter(p,'gpu',gpu_default,@islogical);
            
            parse(p,varargin{:});
        end
        
        function S = ExtractLumpyParam(obj,L_in)
            % Creates a struct with the lumpy background parameters to pass
            % to the mex function
            S = struct;
            S.K   = L_in.K;
            S.dim = uint32(L_in.dim);
            S.B0  = L_in.B0;
            if(length(L_in.b0)~=S.K)
                S.b0  = L_in.b0*ones(S.K,1);
                S.b0 = S.b0(:);
            else 
                S.b0 = L_in.b0(:);
            end
            
            if(all(size(L_in.cov)==[L_in.dim,L_in.dim]))
                % Uniform covariance
                temp = inv(L_in.cov);
                mask = tril(true(size(L_in.cov)));
                temp = temp(mask); 
                cov_vec = temp;
                %unif = 1;
                %nParam = 1+obj.dim*(obj.dim+1)/2;
            elseif(size(L_in.cov,2)==L_in.dim*L_in.K)
                % Non-uniform covariance
                nParam = 1+L_in.dim*(L_in.dim+1)/2;
                cov_vec = zeros(nParam-1,L_in.K);
                mask = tril(true([L_in.dim,L_in.dim]));
                %unif = 0;
                for i=1:obj.K
                    temp = inv(L_in.cov(:,(L_in.dim*(i-1)+1):(L_in.dim*i)));
                    temp = temp(mask);
                    cov_vec(:,i) = temp;
                end
            else
                error('Incorrectly formatted cov! Must be dim-by-dim or dim-by-K*dim.  Size(cov) = (%i,%i)\n dim = %i\n K = %i',size(obj.cov,1),size(obj.cov,2),obj.dim,obj.K);
            end
            
            S.cov = cov_vec(:);
            
%             if(any(size(L_in.cov)~=[S.dim,S.K]))
%                 S.cov = repmat(L_in.cov,[1,S.K]);
%                 S.cov = S.cov(:);
%             else
%                 S.cov = L_in.cov(:);
%             end
            
            S.centers = L_in.centers';  % !!!! Must transpose!!!
            S.centers = S.centers(:);
        end
        
    end
    
end