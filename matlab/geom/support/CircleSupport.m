% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
% 
% File:     CircleSupport.m
% Author:   Nick Henscheid
% Date:     6-2019
% Info:     Class definition file for a circular support set (2D only!!) 
% Inputs:   
% Notes: 
% To Do: 
% 
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

classdef CircleSupport < SupportSet
    
    properties 
       radius = .5; 
       origin = [0,0];
    end
    
    properties (Access = private)
        dim    = 2;  % 2D only!!
    end
    
    methods 
        function obj = CircleSupport(radiusi)
            obj.radius = radiusi;
        end
        
        function plot(obj)
           tplot = linspace(0,2*pi,128);
           X = obj.radius*cos(tplot);
           Y = obj.radius*sin(tplot);
           plot(X,Y,'LineWidth',2.0);
           axis equal;
           xlabel('$x$ (cm)');
           ylabel('$y$ (cm)');
        end
        
        function tau = ExitTimes(obj,r,s)
            % Compute forward and backward exit times
            % r and s are N-by-dim arrays 
            % tau is N-by-2, with tau(:,1)=tau_-,
            % tau(:,2)=tau_+ 
            % if any rays do not intersect, we return NaN
            if(any(size(r)~=size(s)))
                error('r and s must be same size!');
            end
            if(size(r,2)~=2 || size(s,2)~=2)
                error('r and s must be N-by-2');
            end
            N = size(r,1);
            xtilde  = r - repmat(obj.origin,[N,1]);
            xtilde_dot_s = r(:,1).*s(:,1) + r(:,2).*s(:,2);
            norm_xtilde = xtilde(:,1).^2 + xtilde(:,2).^2;
            tau = zeros(N,2);
            tau(:,1) = -xtilde_dot_s + sqrt(xtilde_dot_s.^2 -norm_xtilde + obj.radius^2);
            tau(:,2) = -xtilde_dot_s - sqrt(xtilde_dot_s.^2 -norm_xtilde + obj.radius^2);
            tau(imag(tau)~=0)=nan;
        end
        
        function r = BdyFun(obj,theta)
            % aziumuth angle theta, elevation angle phi
            r = obj.radius;
        end
        
        function a = IsBdy(obj,r)
            a = (norm(r-obj.origin)==obj.radius);
        end
            
        function a = IsBdyPlus(~,r,s)
            a = dot(r,s)>0;
        end
        
        function a = IsBdyMinus(~,r,s)
            a = dot(r,s)<0;
        end
            
        function grid = UnifGrid(~,N)
            % Returns a uniform grid for the support set
            % NEED TO WRITE THIS
            grid = zeros(N);
        end
        
        function samp = UnifSample(~,N)   % Returns N uniformly random samples from the support set
            % NEED TO WRITE THIS
            samp = zeros(N);
        end
        
    end
    
end
