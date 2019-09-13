% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:     RectSupport.m
% Author:   Nick Henscheid
% Date:     6-2019
% Info:     Class definition file for the RectSupport object, a
%           d-dimensional rectalinear support set. 
% Inputs:   
% Notes: 
% To Do:    - allow for arbitrary rotation & origin
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

classdef RectSupport < SupportSet
    properties 
       L = [0,1;0,1]  
       rotation_matrix = eye(2);   
    end
    
    properties (Dependent)
        origin;
        dim; 
    end
    
    methods 
        function obj = RectSupport(varargin)
            p = inputParser;
            default_L = [0,1;0,1];
            isvalidL  = @(x) size(x,2)==2;
            addParameter(p,'L',default_L,isvalidL);
            p.parse(varargin{:});
            obj.L = p.Results.L;
        end
        
        function val = get.dim(obj)
            val = size(obj.L,1);
        end
        
        function val = get.origin(obj)
            val = zeros(1,size(obj.L,1));
            for i=1:size(obj.L,1)
                val(i) = obj.L(i,1);
            end
        end
              
        function plot(obj,varargin)
            % alpha is the opacity 
            p = inputParser; 
            addParameter(p,'axes',[]);
            addParameter(p,'alpha',0.5);
            p.parse(varargin{:});
            
            alpha = p.Results.alpha;
            
            l = obj.L;
            if(obj.dim==2)
               x = [l(1,1),l(1,2),l(1,2),l(1,1),l(1,1)];
               y = [l(2,1),l(2,1),l(2,2),l(2,2),l(2,1)];
               line(x,y,'LineWidth',2);
               xlabel('$x$ (cm)');ylabel('$z$ (cm)');
            elseif(obj.dim==3)
               lx = l(1,2)-l(1,1);
               ly = l(2,2)-l(2,1);
               lz = l(3,2)-l(3,1);
               plotcube([lx,ly,lz],obj.origin,alpha,[0.1,0.1,0.2]);
               xlabel('$x$ (cm)');ylabel('$y$ (cm)');zlabel('$z$ (cm)');  
            end
            axis equal;
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
            tau = 0;
        end
      
        
        function X = UnifGrid(obj,N,varargin)
            p = inputParser;
            addParameter(p,'centered',0);
            parse(p,varargin{:});
            l = obj.L;
            if(p.Results.centered)
                dx = (l(1,2)-l(1,1))/N(1);
                dy = (l(2,2)-l(2,1))/N(2);
                xtmp = linspace(l(1,1)+dx/2,l(1,2)-dx/2,N(1));
                ytmp = linspace(l(2,1)+dy/2,l(2,2)-dy/2,N(2));
            else
                xtmp = linspace(l(1,1),l(1,2),N(1));
                ytmp = linspace(l(2,1),l(2,2),N(2));
            end
            
            if(obj.dim ==2)
                [xx,yy] = meshgrid(xtmp,ytmp);
                X = cat(3,xx,yy);
            elseif(obj.dim == 3)
                if(p.Results.centered)
                    dz = (l(3,2)-l(3,1))/N(3);
                    ztmp = linspace(l(3,1)+dz/2,l(3,2)-dz/2,N(3));
                else
                    ztmp = linspace(l(3,1),l(3,2),N(3));
                end
                [xx,yy,zz] = meshgrid(xtmp,ytmp,ztmp);
                X = cat(4,xx,yy,zz);
            end
        end
        
        function X = UnifSample(obj,N)
            X = rand(N,obj.dim);
            l = obj.L; 
            dx = l(1,2)-l(1,1);
            dy = l(2,2)-l(2,1);
            X(:,1) = obj.L(1,1) + dx*X(:,1);
            X(:,2) = obj.L(2,1) + dy*X(:,2);
            if(obj.dim==3)
                dz = l(3,2)-l(3,1);
                X(:,3) = obj.L(3,1) + dz*X(:,3);
            end
        end    
    end
end
