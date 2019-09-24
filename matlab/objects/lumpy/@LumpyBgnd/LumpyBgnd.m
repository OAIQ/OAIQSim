% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:     LumpyBgnd.m
% Author:   Nick Henscheid
% Date:     9-2016, 9-2019
% Info:     u = LumpyBgnd(varargin). The lumpy background object.
%           The name-value pairs in varargin determine the random field     
% Inputs: 
%           'b0' (a scalar)   The lump amplitude l(x) = b0*l0(x)
%           'B0' (a scalar)   The DC offset i.e. f(x) = B0 + sum()
%           'cov' (a scalar or dim-by-1 vector)
%           'Kbar' (a positive scalar) 
%           'gpu',(0 or 1) If gpu = 1 and a gpu is available, the CUDA 
%                          version is used.
%           'N' (an integer or dim-by-1 vector of integers)
% Notes:   
% To Do: 
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

classdef LumpyBgnd < handle
    properties (SetObservable = true)
        b               % Lump "Amplitude"
        B0              % DC offset
        cov             % Lump covariance matrix 
        Kbar            % Mean number of lumps
        centers         % Lump centers
        gpu             % Evaluate using gpu or not
        N               % Default number of eval. pts in ea. direction 
        K_distr;        % Probability distribution for the # of lumps
        centers_distr;  % Probability distribution for the lump centers 
        b_distr;        % Probability distribution for the lump amplitudes
    end
    
    properties (Dependent)
        L;    % Bounding box for evaluation (lump centers can extend slightly beyond to avoid edge effect issues)
        dim; 
        K; 
    end

    properties (SetAccess=private)
        support;    % Support set for evaluation.  Can only set when LumpyBgnd is first initialized!
        pad_factor = 3;      % Padding factor (so that boundary effects don't occur)
        show_warnings = 1;
        random_number_of_lumps = 1; % Randomize number of lumps when randomize() is called?
        normalize;                  % Produce an approximate PDF?      
    end
    
    % Standard methods
    methods
        function obj = LumpyBgnd(varargin)
            p = obj.ParseInputs(varargin{:});
            obj.Kbar     = p.Results.Kbar;
            %obj.b        = p.Results.b;
            obj.B0       = p.Results.B0;
            obj.support  = p.Results.support;
            if(numel(p.Results.cov)==1)
                obj.cov = p.Results.cov*eye(obj.support.dim);
            else
                obj.cov = p.Results.cov;
            end
            obj.centers = p.Results.centers;
            obj.K_distr = p.Results.K_distr;
            obj.centers_distr = p.Results.centers_distr;
            obj.b_distr = p.Results.b_distr;
            
            if(gpuDeviceCount()>0)
                obj.gpu     = p.Results.gpu;
            else
                obj.gpu     = 0;
            end
            obj.N       = p.Results.N;
            if(numel(obj.centers)==0)
                % Need to generate random lump centers
                obj.randomize();
            end
            addlistener(obj,'Kbar','PostSet',@LumpyBgnd.handlePropertyEvents);
            addlistener(obj,'cov','PostSet',@LumpyBgnd.handlePropertyEvents);
        end
        
        % Get and set methods for dependent properties
        function val = get.L(obj)
            val = obj.support.L;
        end
        
        function val = get.dim(obj)
            val = obj.support.dim;
        end
        
        function val = get.K(obj)
            val = size(obj.centers,1);
        end
        
        function set.K(obj,value)
            if(isnumeric(value)&&value>=0&&(mod(value,1)==0))
                obj.K = value;
            else
                error('Invalid K value');
            end
        end
        
        function TurnOffWarnings(obj)
            obj.show_warnings = 0;
        end
        
        function TurnOnWarnings(obj)
            obj.show_warnings = 1;
        end
        
        function SetPadFactor(obj,x)
            obj.pad_factor = x;
        end
        
        function SetRandomLumps(obj,x)
            obj.random_number_of_lumps = x;
        end
        
        function set.gpu(obj,value)
            if(isnumeric(value)||islogical(value))
                if(value==1||value==true)
                    if(gpuDeviceCount>0)
                        obj.gpu = value;
                    end
                elseif(value==0||value==false)
                    obj.gpu = value;
                else
                    error('Invalid value for gpu!');
                end
            else
                error('Invalid value for gpu!');
            end
        end
        
        % Externally defined functions
        p = ParseInputs(varargin);
        u = eval(obj,X,XSize);
        obj = randomize(obj,varargin)
        varargout = plot(obj,varargin);
        z = minus(x,y);  % Method to subtract two lumpy backgrounds
        z = plus(x,y);   % Method to add two lumpy backgrounds
        
        function U = sample(obj,Ns)
            % Generates Ns samples of the lumpy background, returning it in
            % a cell array U (U{i} is a sample for each i)
            U = cell(Ns,1);
            for i=1:Ns
                progressbar(i/Ns);
                obj.randomize;
                U{i} = obj.eval;
            end
        end
    end
    
    % Static methods
    methods (Static)
        handlePropertyEvents(src,evnt);   % Defined externally
    end
    
end