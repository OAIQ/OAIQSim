% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File: PPDetector.m
% Author:  Nick Henscheid
% Date:    9-2017, 5-2018, 6-2019
% Info:    Matlab class for planar photon-processing scintillation
%          detectors
% Inputs:  
% Notes:   - This is a model for a planar scintillation detector that
%            estimates photon interaction position (x) at the face of the
%            detector
%          - Detector is modeled as a slab with dimensions [w,d] (2d)
%            or [w,h,d] (3d)
%          - Collimator is modeled as a perfect hard-edge pinhole of width
%            w, assumed centered on the axis of the detector at distance p
%            from the detector face.
%          - Estimation blur is assumed to be unbiased and Gaussian with
%            covariance matrix C
%
% To Do: 
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------
classdef PPDetector
    
    properties 
        detl;      % Dimensions of detector; detl = [w,d] in 2d, 
                   % detl = [w,h,d] in 3d.
        detpos;    % Detector origin relative to system coordinate origin
                   % if 2d: detpos = [rho,theta]  (standard polar rep.)
                   % if 3d :
                   % detpos = [rho,phi,z], [x,y,z], [rho,phi,theta]
                   % for resp. cyl., cart., and sph. coords
                   % See doc for description.
        npin;      % Number of pinholes
        pinw;      % Pinhole widths  (vector size npin-x-1)
        pinc;      % Pinhole centers (vector size npin-x-1)
        pind;      % Pinhole distance to detector face
        %mdrf;     % Not used yet. 
        blursigma; % 
        dim;       % Dimension (2 or 3)
        coordsys;  % Coordinate system ('cart', 'cyl' or 'sph')
    end

    methods
        function obj = PPDetector(varargin)
            p = obj.ParseInputs(varargin{:});
            obj.dim        = p.Results.dim;
            obj.detl       = p.Results.detl;
            obj.detpos     = p.Results.detpos;
            obj.npin       = length(p.Results.pinc);
            obj.pinw       = p.Results.pinw;
            obj.pinc       = p.Results.pinc;
            obj.pind       = p.Results.pind;
            obj.blursigma  = p.Results.blur;
            obj.coordsys   = p.Results.coordsys;
        end
        %%% Externally defined functions %%%
        plot(obj); % Detector plotting function
        p = ParseInputs(obj,varargin); % Input parser
        beta = getVisibleLimits(obj,bdyfun); % Get visible limits given a support object
        beta = getPinholeAngles(obj,x); % Get pinhole angles
        plotFOV(obj,x,support); % Plot field of view for the detector   
        M = SimulateMDRF(obj);   % MDRF Simulator 
    end
    
end