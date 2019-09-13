% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:     BinnedModeData.m
% Author:   Nick Henscheid
% Date:     6-2019
% Info:     Class definition file for binned-mode data.  Binned data
%           assumes a real or artificial pixelation/binning of the
%           detector.  
% Inputs:   
% Notes:    
% To Do: 
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------
classdef BinnedModeData
    
    properties 
        num_detetector;  % Number of detectors
        pixel_size;      % [Lx,Ly] 
        num_pixels;      % [px,py]
        g;               % data array
    end
    
    
    
    methods 
        function obj = BinnedModeData
            
        end
        
    end
end