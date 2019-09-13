% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:    FASTSPECTSystem.m
% Author:  Nick Henscheid
% Date:    9-2016, 6-2019
% Info: A matlab class for FASTSPECT-type fixed detector system simulation.  
% Inputs:  - 
%          - 
% Notes: 
% To Do: 
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------
classdef FastSPECTSystem < handle
    properties 
        dets;        % Array of detector objects
        det_type;     % 'pp' or 'pix'
        support;     % Support object   
        %supporttype; % 'sph' = Sphere, 'cart' = rectalinear
                     % or 'cyl' = Cylindrical object support.  
                     % Don't think this is necessary!
    end
    
    properties (Dependent)
        dim;         % Dimension (2d or 3d) 
        num_det;        % Number of detectors
    end
    
    properties (SetAccess = private)
        transport_solver;  % The RTE solver 
        detector_solver;   % The detector simulator
    end
    
    methods
        function obj = FastSPECTSystem(varargin)
            p            = obj.ParseInputs(varargin{:});
            obj.dets     = p.Results.dets;
            obj.support  = p.Results.support;
            obj.det_type = p.Results.det_type;
        end
        
        % Get and set methods for dependent properties 
        function val = get.dim(obj)
            val = obj.support.dim;
        end
        
        function val = get.num_det(obj)
            val = length(obj.dets);
        end
        
        % Externally defined functions 
        plotMeasurementSet(obj,gamma,j);
        p    = ParseInputs(varargin);
        gbar = SimulateMeanData(obj,f);   % System simulator (mean)
        g    = SimulateRandomData(obj,f); % System simulator (random)
        M    = SimulatePSF(obj);          % PSF simulator 
        
        % Functions to modify detector properties
        function SetPinholeWidth(obj,w)
            % Sets a common pinhole size for all the detectors
            for i=1:obj.num_det
                obj.dets(i).pinw = w;
            end
        end
        
        function SetPinholeDist(obj,d)
            % Sets a common pinhole distance for all the dets
            for i=1:obj.num_det
                obj.dets(i).pind = d;
            end
        end
        
        function SetDetSize(obj,L)
            % Sets a common detector size 
            for i=1:obj.num_det
                obj.dets(i).detl = L;
            end
        end
        
        % Plotting functions
        function plot(obj)
            obj.plotSupport;
            for i=1:obj.num_det
                obj.dets(i).plot;
                hold on;
            end
            hold off;
            title('System Geometry');
        end
      
        function plotSupport(obj)
            obj.support.plot;
        end
        
        function plotFOV(obj,x,j)
            if nargin<3
                j = 1; 
            end
            if nargin<2
                x = 0;
            end
            figure;
            plot(obj);
            for i=1:length(j)
                plotFOV(obj.dets(j(i)),x,obj.support)
            end
            title(sprintf('Adjoint Ray Bundle for %1.3f$\\leq x\\leq$ %1.3f',min(x),max(x)),'interpreter','latex','FontSize',14);
        end
        
        function ans = isAngleBetween(obj,beta1,beta2,beta)
            % Checks if an angle beta is between two angles beta1 and beta2
            % mod 2*pi
            beta = beta-beta1; 
            if(beta<0)
                beta = beta+2*pi;
            end
            beta2 = beta2-beta1;
            if(beta2<0)
                beta2 = beta2+2*pi;
            end
            if(beta<=beta2)
                ans = 1;
            else
                ans = 0;
            end
        end
        
        function plotVisibleBdy(obj,j)
            % For each detector, find the portion of the boundary that is
            % accessible and plot it 
            if nargin<2
                j = 1:obj.num_det;
            end
            figure;
            plot(obj);hold on;
            for i=1:length(j)
                beta = obj.dets(j(i)).getVisibleLimits(@(x)obj.support.bdyFun(x));  % Beta will be 2-by-npin
                for k=1:obj.dets(j(i)).npin
                    tplot = linspace(min(beta(:,k)),max(beta(:,k)),100);
                    rplot=obj.support.bdyFun(tplot);
                    plot(rplot.*cos(tplot),rplot.*sin(tplot),'r','LineWidth',2.5);
                end
            end
            hold off;
        end
    end
end