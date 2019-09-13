% ---------------------------------------------------------------------------
% **********************CGRI SIMULATION TOOLBOX***************************
% u = LumpyDiffusion(). The lumpy background object.

%  **** UPDATE!!!! ****
%  **** UPDATE!!!! ****
%  **** UPDATE!!!! ****
% u = LumpyBgnd(X,varargin) evaluates a lumpy background image at the
% coordinates in the N-by-d matrix X.
% The name-value pairs in varargin determine if the background is generated
% from scratch or using pre-determined lump centers.
%   ex: [u,c] = LumpyBgnd(X,'Kbar',100,'b0',1,'C0',0,'rb2',0.001); 
%   generates a random lumpy background from scratch with the parameters 
%   (Kbar,b0,C0,rb2) = (100,1,0,0.001), evaluted at the coordinates X.  X
%   is N-by-d, with each row being an evaluation point (x,y) or (x,y,z)
%   ex: u = LumpyBgnd(X,'centers',c,'b0',1,'C0',0,'rb2',0.001); evalutes a
%   lumpy background with pre-determined lump centers c (K-by-d) and lump
%   parameters (b0,C0,rb2);
% If gpu = 1 and a gpu is available, the CUDA version is used.
% 

% File: LumpyDiffusion.m
% Author:  Nick Henscheid
% Date:    8-2017
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without technical
% support, and with no warranty, express or implied, as to its usefulness for
% any purpose.
% ---------------------------------------------------------------------------

classdef LumpyDiffusion < handle
    properties (SetObservable = true)
        b0      % Lump "Amplitude" vector
        C0      % DC offset
        rb2     % Lump "width" vector
        Kbar    % Mean number of lumps
        centers % Emission points
        times   % Emission times
        gpu     % Evaluate using gpu or not
        nx      % number of pixels in ea. direction (assuming square/cube)
        pdf     % forces a normalization where the resulting samples are probability densities
        D       % Diffusion coefficient
        t       % Evaluation time
        auc     % Evaluate AUC(T) instead?
        K       % Actual number of lumps
        lambda  % Intensity function - must be a function handle @(x) lambda(x)
                %                      x(:,1:dim) are positions, x(:,dim+1)
                %                      is time
        tmax    % Max time for samples 
    end
   
    properties (SetAccess = private)
        dim     % Dimension
    end
    
    methods
        function obj = LumpyDiffusion(varargin)
            p = obj.parseinputs(varargin{:});
            obj.dim     = p.Results.dim;
            obj.Kbar    = p.Results.Kbar;
            obj.b0      = p.Results.b0;
            obj.C0      = p.Results.C0;
            obj.rb2     = p.Results.rb2;
            obj.centers = p.Results.centers;
            obj.times   = p.Results.times;
            obj.D       = p.Results.D; % Diffusion coef
            if(sign(gpuDeviceCount))
                obj.gpu     = p.Results.gpu;
            else
                obj.gpu     = 0;
            end
            obj.nx      = p.Results.nx;
            obj.pdf     = p.Results.pdf;  % Normalize so the result is a kernel density estimator
            obj.t       = p.Results.t;
            obj.auc     = p.Results.auc;
            obj.lambda  = p.Results.lambda;
            obj.tmax    = p.Results.tmax;
            if(any(isnan(obj.centers(:)))||any(isnan(obj.times(:))))
                % Need to generate centers and times 
                obj.randomize();
            end
            %obj.K = size(obj.centers,1);
            if(obj.pdf)
                obj.rb2 = 2*obj.rb2;
                obj.b0 = sqrt(pi*obj.rb2);
                obj.C0 = 0;
            end
        end
        
        function obj = randomize(obj)
            % Want to convert this to using PPP to generate points using
            % the intensity function lambda. 
            
            if(obj.dim==2)
                ntemp = 128;
                [xtemp,ytemp] = meshgrid(linspace(0,1,ntemp));
                X = [xtemp(:),ytemp(:),linspace(0,obj.tmax,ntemp^2)'];
                lambda_0 = obj.lambda(X); upper_bound = 1.2*max(lambda_0(:));
                pad = 2*sqrt(obj.D*obj.tmax);
                Z =  sortrows(PPP(obj.lambda,obj.Kbar,3,[-pad,1+pad;-pad,1+pad;0,obj.tmax],upper_bound,50),3);
            elseif(obj.dim==3)
                ntemp = 128;
                [xtemp,ytemp,ztemp] = meshgrid(linspace(0,1,ntemp));
                X = [xtemp(:),ytemp(:),ztemp(:),linspace(0,obj.tmax,ntemp^3)'];
                lambda_0 = obj.lambda(X); upper_bound = 1.2*max(lambda_0(:));
                pad = 2*sqrt(obj.D*obj.tmax);
                Z =  sortrows(PPP(obj.lambda,obj.Kbar,4,[-pad,1+pad;-pad,1+pad;-pad,1+pad;0,obj.tmax],upper_bound,50),4);
            end
            obj.K = length(Z);
            obj.times   = Z(:,obj.dim+1);
            obj.centers = Z(:,1:obj.dim);
            
            
%             rng('shuffle');
%             pad = 3*sqrt(obj.rb2/2);  %3-sigma padding to improve stationarity
%             minx = -pad;  %Note: assuming supp(u) \subset [0,1]^dim
%             maxx = 1+pad;
%             miny = -pad;
%             maxy = 1+pad;
%             if(obj.dim==3)
%                 minz = -pad;
%                 maxz = 1+pad;
%             end
%             obj.K = poissrnd(obj.Kbar);
%             obj.times = rand(obj.K,1);
%             obj.centers = zeros(obj.K,obj.dim);
%             obj.centers(:,1) = minx + (maxx-minx).*rand(obj.K,1);
%             obj.centers(:,2) = miny + (maxy-miny).*rand(obj.K,1);
%             if(obj.dim==3)
%                 obj.centers(:,3) = minz + (maxz-minz).*rand(obj.K,1);
%             end
            
            disp(size(obj.times)); disp(obj.K);disp(size(obj.centers));
        end
        
        function u = eval(obj,X,XSize)
            % This function will evaluate the texture for the sample points
            % defined in the array X at time t>0.
            % This should be updated once the GridData object is ready to
            % go!
            if(obj.K ~= size(obj.centers,1))
                error('Lump number and size of centers are not equal!');
            end
            if(size(obj.centers,1)~=size(obj.times,1))
                error('Number of centers and number of times are not equal!');
            end
            
            if(nargin<2) % No array provided
                if(obj.dim==2)
                    [xtemp,ytemp] = meshgrid(linspace(0,1,obj.nx));
                    X = [xtemp(:),ytemp(:)];
                    XSize = [obj.nx,obj.nx];
                elseif(obj.dim==3)
                    [xtemp,ytemp,ztemp] = meshgrid(linspace(0,1,obj.nx));
                    X = [xtemp(:),ytemp(:),ztemp(:)];
                    XSize = [obj.nx,obj.nx,obj.nx];
                end
            end
            X = X';
            nEval = size(X,2)
            c_vec = obj.centers';
            t_vec = obj.times';
            if(obj.gpu)
                disp('Calculating with GPU...')
                tic
                u = lumpy_diff_mex_gpu(uint32(obj.dim),uint32(obj.K),uint32(nEval),uint32(obj.auc),single(c_vec(:)),single(X(:)),single(t_vec(:)),single(obj.t),single(obj.C0),single(obj.D));
                disp(['gpu time = ',num2str(toc)])
                disp(['length(u) = ',num2str(length(u))])
            else
                disp('Calculating with CPU...')
                u = lumpy_diff_mex_cpu(uint32(obj.dim),uint32(obj.K),uint32(nEval),uint32(obj.auc),single(c_vec(:)),single(X(:)),single(t_vec(:)),single(obj.t),single(obj.C0),single(obj.D));
                disp('Done!')
            end
            if(nargin<2)  % No array provided, reshaping to locally generated grid size
                u = reshape(u,XSize);
            end
            if(nargin==3) % User supplied a grid size, reshape to it
                u = reshape(u,XSize);
            end
            if(obj.pdf)
                disp('pdf!');
                u = u/obj.K;
            end
        end
        
        function k = corrfun(obj,dr)
            % Computes the correlation function at a given delta-r 
            k = (obj.Kbar*obj.b0^2/(2*pi*obj.rb2^2))*exp(-dr.^2/(2*obj.rb2^2)); 
        end
        
        function U = sample(obj,Ns)
            % Generates Ns samples of the lumpy background, returning it in
            % a cell array U (U{i} is a sample for each i)
            U = cell(Ns,1);
            for i=1:Ns
                obj.randomize;
                U{i} = obj.eval;
            end
        end
                
        function p = parseinputs(~,varargin)
            p = inputParser;
            p.CaseSensitive = true;
            Kbar_default    = 100;
            b0_default      = 1;
            C0_default      = 0;
            rb2_default     = 0.005;
            dim_default     = 2;
            centers_default = NaN;  %??  Why?
            times_default   = NaN;  %?? 
            gpu_default     = 0;
            nx_default      = 64;
            auc_default     = 0;
            D_default       = 1e-4;  % For molecular diffusion
            t_default       = 1;  % Evaluation time
            tmax_default    = 1;
            lambda_default  = @(x) 1;   % Mean function
            isnonneg        = @(x) (x>=0);
            ispos           = @(x) (x>0);
            isvaliddim      = @(x) (x==2||x==3);
            isbinary        = @(x) (x==0||x==1);
            ishandle        = @(x) strcmp(class(x),'function_handle');
            addParameter(p,'b0',b0_default,ispos);
            addParameter(p,'C0',C0_default,isnonneg);
            addParameter(p,'rb2',rb2_default,ispos);
            addParameter(p,'dim',dim_default,isvaliddim);
            addParameter(p,'Kbar',Kbar_default,ispos);
            addParameter(p,'centers',centers_default,@isnumeric);
            addParameter(p,'times',times_default,@isnumeric);
            addParameter(p,'gpu',gpu_default,isbinary);
            addParameter(p,'nx',nx_default,ispos);
            addParameter(p,'pdf',false,@islogical);
            addParameter(p,'D',D_default,ispos);
            addParameter(p,'t',t_default,ispos);
            addParameter(p,'auc',auc_default,isbinary);
            addParameter(p,'lambda',lambda_default,ishandle);
            addParameter(p,'tmax',tmax_default,ispos);
            parse(p,varargin{:});
        end
    end  
end