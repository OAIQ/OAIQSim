% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
% 
% File:     CylinderSupport.m
% Author:   Nick Henscheid
% Date:     6-2019
% Info:     Class definition file for cylindrical support set (3D only!!)
% Inputs:   
% Notes: 
% To Do:    - allow arbitrary axis
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------


classdef CylinderSupport < SupportSet
    % Cylindrical support set with axis aligned to z
    % 
    properties 
       radius = .5; 
       length = 1;
       origin = [0,0,0];
       axis   = [0,0,1];  % Not currently used
    end
    
    properties (Access = private)
       dim    = 3; % 3D only!!
    end
    
    methods 
        function obj = CylinderSupport(radiusi,lengthi)
            obj.radius = radiusi;
            obj.length = lengthi;
        end
        
        function plot(obj)
           [X,Y,Z] =  cylinder(1,128);   
           X = obj.radius*X + obj.origin(1);          
           Y = obj.radius*Y + obj.origin(3);
           Z = obj.length*Z + obj.origin(2);
           surface(X,Z,Y);
           axis equal;
           xlabel('$x$ (cm)');
           ylabel('$z$ (cm)');
           zlabel('$y$ (cm)');
           % Note: Z and Y are flipped so that Z-axis is the patient bed
        end
        
        function r = BdyFun(obj,theta,z)
            % aziumuth angle theta, elevation angle phi
            r = obj.origin + [obj.radius*cos(theta),obj.radius*sin(theta),z];
        end
        
        function a = IsBdy(obj,r)
            a = (norm(r-obj.origin)==obj.radius);
        end
            
        function a = IsBdyPlus(obj,r,s)
            a = dot(r,s)>0;
        end
        
        function a = IsBdyMinus(obj,r,s)
            a = dot(r,s)<0;
        end
            
        function tau = ExitTimes(obj,r,s)
            error('Not implemented yet!');
            tau = NaN;
        end
        
        function grid = UnifGrid(~,N)
            % Returns a uniform grid for the support set
            grid = zeros(N);
        end
        
        function samp = UnifSample(~,N)   % Returns N uniformly random samples from the support set
            samp = zeros(N);
        end
    end
    
end
