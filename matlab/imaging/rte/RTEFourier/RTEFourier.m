% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:     RTEFourier.m
% Author:   Nick Henscheid
% Date:     9-2016, 6-2019
% Info:     H = RTEFourier(varargin).  A ballistic RTE solver for Fourier 
%           modes.
% Inputs:   The RTEFourier object has the following optional inputs: 
%           - 'k' is the mode vector i.e. k = [k1,k2] (2D) or k =
%           [k1,k2,k3] (3D)
%           - 'cm' is the medium speed of light
%           - 'hmin' is the ballistic RTE integrator minimum step size
% Notes:  
% To Do:    - Scattering
% 
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------
classdef RTEFourier < handle
    properties (SetObservable = true)
       k;     % Fourier mode vector
       cm;    % Light speed
       hmin;  % Minimum step length for RTE solver
    end
    
    properties (SetAccess = private)
       gpu;
    end
    
    properties (Dependent)
        dim;
    end
    
    methods
        function obj = RTEFourier(varargin)
            p         = obj.parseinputs(varargin{:});
            obj.k     = p.Results.k(:)';  
            obj.cm    = p.Results.cm;
            obj.hmin  = min(p.Results.hmin,1/(2*norm(obj.k)));
            hasgpu = logical(gpuDeviceCount);
            if(hasgpu)
                obj.gpu   = p.Results.gpu;
            else
                disp('No GPU found.  Using CPU!');
                obj.gpu   = false;
            end
        end
        
        
        function val = get.dim(obj)
            val = numel(obj.k);
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
        
        function w = eval(obj,R,S,E)
           % Evaluates w for vectors R,S and E.  R,S and E must be the same
           % length, dim*neval-by-1.  length(w) = neval.
           if(isequal(size(R),size(S)))
               if(size(R,2)>1||size(S,2)>1||size(E,2)>1)
                   error('Wrong size inputs! R,S and E must be the same length, dim*neval-by-1.  length(w) = neval')
               else 
                   neval = length(R)/obj.dim;
                   if(obj.gpu)
                       disp('Computing with GPU');
                       % Should compile two versions!
                       w = rte_fourier_mex_gpu(uint32(obj.dim),R,S,obj.k,uint32(neval),obj.cm,obj.hmin);
                       %w = rte_fourier_mex_gpu(uint32(obj.dim),single(R),single(S),single(obj.k),uint32(neval),single(obj.cm),single(obj.hmin));
                   else
                       disp('Computing with CPU');
                       disp(obj)
                       nstep = uint32(ceil(1/obj.hmin));
                       w = rte_fourier_mex_cpu(uint32(obj.dim),R,S,obj.k,uint32(neval),obj.cm,nstep);
                   end
               end
           else
               error('Inputs must be same size');
           end
        end

        function p = parseinputs(~,varargin)
            p = inputParser;
            p.CaseSensitive = true;
            k_default       = [0,0];
            cm_default      = 299792458;
            hmin_default   = 0.001;   % Step size for integrator
            gpu_default     = false;
            kcheck = @(k)(isnumeric(k)&&(numel(k)==2||numel(k)==3));
            
            addParameter(p,'k',k_default,kcheck);
            addParameter(p,'cm',cm_default,@isnumeric);
            addParameter(p,'hmin',hmin_default,@isnumeric);
            addParameter(p,'gpu',gpu_default,@islogical);
            
            parse(p,varargin{:});
        end
    end
    
end