% ---------------------------------------------------------------------------
% **********************CGRI SIMULATION TOOLBOX***************************
% RTEFourier. A ballistic RTE solver for Fourier modes.

%  **** UPDATE ****
% T = RTEFourier

% File: LumpyBgnd.m
% Author:  Nick Henscheid
% Date:    9-2016
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without technical
% support, and with no warranty, express or implied, as to its usefulness for
% any purpose.
% ---------------------------------------------------------------------------
classdef RTEVoxel < handle
    properties (SetObservable = true)
       k;     % Fourier mode vector
       cm;    % Light speed
       hmin; % Number of RK4 steps
    end
    
    properties (SetAccess = private)
       gpu;
    end
    
    methods
        function obj = RTEFourier(varargin)
            p         = obj.parseinputs(varargin{:});
            obj.k     = p.Results.k;
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
           % length, neval-by-1.  length(w) = neval.
           if(isequal(size(R),size(S),size(E)))
               if(size(R,2)>1||size(S,2)>1||size(E,2)>1)
                   error('Wrong size inputs!')
               else 
                   neval = length(R)/2;
                   if(obj.gpu)
                       disp('Computing with GPU');
                       w = rte_fourier_mex_gpu(R,S,obj.k,uint32(neval),obj.cm,obj.hmin);
                   else
                       disp('Computing with CPU');
                       disp(obj)
                       nstep = uint32(ceil(1/obj.hmin))
                       w = rte_fourier_mex_cpu(R,S,obj.k,uint32(neval),obj.cm,nstep);
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
            hmin_default   = 0.01;   % Step size for integrator
            gpu_default     = false;
            
            addParameter(p,'k',k_default,@isnumeric);
            addParameter(p,'cm',cm_default,@isnumeric);
            addParameter(p,'hmin',hmin_default,@isnumeric);
            addParameter(p,'gpu',gpu_default,@islogical);
            
            parse(p,varargin{:});
        end
    end
    
end