% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:    VoxelData.m
% Author:  Nick Henscheid
% Date:    10/2017, 5/2018, 8/2018, 6/2019
% Info:    Class definition file for VoxelData, a class for
%          gridded/voxelized data.  
% Inputs: 
% Notes:   - More appropriate for voxelized *objects*; for pixelated/binned
%            image data, use the BinnedModeData object.
%          - Assumes a uniform rectilinear voxel grid
%          - the size of the voxels is set by the grid limits and the
%          number of voxels in each dimension i.e. if the grid is on the
%          unit square (in cm) and [n,m] = [128,64], the voxels are
%          1e-2/128 m = 78.25 microns by 1e-2/64 = 156.25 microns.
%          - The voxel values can either represent function samples f(x_j),
%          with x_j typically at the voxel center, or represent voxel 
%          averages (integrals of f(x) over V_j) or probabilities.
% To Do:   - GUI 
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------
classdef VoxelData < handle      
    properties  (SetObservable = true)
       support RectSupport; % Rectangular support set 
       data_array;          % Voxel data values, arranged as n-by-m-by-datadim or 
                            % l-by-n-by-m-by-datadim 
    end
    
    properties (SetAccess = private)
       imfig;    % To store a display figure (MVC DESIGN...????) 
    end
    
    properties (Dependent)
        grid_size;        % Grid size, e.g. [Nx,Np], [Nx,Ny,Np] or [Nx,Ny,Nz,Np]
        dim;      % Grid dimension, 1-3 .  griddim = length(gridsize);
    end
    
    methods
        function obj = VoxelData(varargin)
            inputs      = obj.ParseInputs(varargin{:});    
            obj.support = inputs.support; 
            obj.data_array       = double(inputs.data_array);
        end

        function val = get.grid_size(obj)
            sz = size(obj.data_array);
            if(length(sz)>obj.dim)
                val = sz;
            elseif(length(sz)==obj.dim)
                val = [sz,1];
            else
                error('Something wrong with dimensions!')
            end
        end
        
        function val = get.dim(obj)
            val = obj.support.dim;
        end
        
        function z = unroll(obj)
            % Currently only works for p = 1...
            if(obj.dim==2)
                xx = obj.X(:,:,1);
                yy = obj.X(:,:,2);
                z = [xx(:),yy(:)];
            end
            if(obj.dim==3)
                xx = obj.X(:,:,:,1);
                yy = obj.X(:,:,:,2);
                zz = obj.X(:,:,:,3);
                z  = [xx(:),yy(:),zz(:)];
            end
        end
        
        % Declare input parser
        p = ParseInputs(~,varargin);
        
        % Declare slice method
        %slice(obj);
        
        % Imagesc method
        imagesc(obj);
        
        % Mean method
        function y = mean(obj)
            y = mean(obj.data_array(:));
        end
        
        % A function to plot the voxel grid
        function PlotGrid(obj,xx,yy)
            figure(obj.imfig);
            hold on;           
            for ii=1:obj.grid_size(1)
                line([xx(ii),xx(ii)],[obj.support.L(2,1),obj.support.L(2,2)],'Color','k','LineWidth',0.25);
            end
            
            for jj=1:obj.grid_size(2)
                line([obj.support.L(1,1),obj.support.L(1,2)],[yy(jj),yy(jj)],'Color','k','LineWidth',0.25);
            end
            hold off;
        end 
        
        function plot(obj,varargin)
            p = inputParser;
            addParameter(p,'grid',0);
            addParameter(p,'axes',[]);
            addParameter(p,'alpha',0.5);
            parse(p,varargin{:});
            
            obj.support.plot('axes',p.Results.axes,'alpha',p.Results.alpha); 
            obj.imfig = gcf;
            hx = (obj.support.L(1,2)-obj.support.L(1,1))/(obj.grid_size(1));
            xx = (obj.support.L(1,1)+hx/2):hx:(obj.support.L(1,2)-hx/2);
            hy = (obj.support.L(2,2)-obj.support.L(2,1))/(obj.grid_size(2));
            yy = (obj.support.L(2,1)+hy/2):hy:(obj.support.L(2,2)-hy/2);
            hold on;
            imagesc(xx,yy,obj.data_array);
            axis image;
            hold off;
            if(p.Results.grid)
                xx = xx-hx/2;xx(end+1) = obj.support.L(1,2);
                yy = yy-hy/2;yy(end+1) = obj.support.L(2,2);
                obj.PlotGrid(xx,yy); 
            end
        end
        
        function gridCallback(obj,src,evt)
            chillens = allchild(obj.imfig.Children(end));
            if(length(chillens)==1)
                obj.plotGrid;
            end
            if(src.Value==0)
                set(chillens(1:end-1),'Visible','off')
            else
                set(chillens(1:end-1),'Visible','on')
            end        
        end
         
        function paramCallback(obj,src,evt)
            selection = src.Value;
            obj.imfig.Children(end).Children(end).CData = obj.data_array(:,:,selection);
        end
      
        function imfigCloseFcn(obj,src,data)
            delete(obj.imfig);
            obj.imfig = [];
        end
        
        % Declare slice callback
        slicecallback(obj,source,event,panel,handles);
        
        function isocallback(obj,source,event,panel,handles)
            isolev = handles.isolev.Value;set(handles.isotxt,'String',['isolevel= ',num2str(isolev)]);
            ax = gca;
            delete(ax.Children(:))
            %delete(ax.Children(3));
            fprintf('Rendering isosurface...')
            hpatch = patch(isosurface(obj.X(:,:,:,1),obj.X(:,:,:,2),obj.X(:,:,:,3),obj.data_array,isolev));
            isonormals(obj.X(:,:,:,1),obj.X(:,:,:,2),obj.X(:,:,:,3),obj.data_array,hpatch)
            hpatch.FaceColor = 'k';
            hpatch.FaceAlpha = 0.15;
            hpatch.EdgeColor = 'none';
            fprintf('done!\n')
            obj.slicecallback(source,event,panel,handles);
        end 
    end
    
    methods (Access = protected)
       
    end
    
    
end
    