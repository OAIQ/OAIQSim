% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:     SupportSet.m
% Author:   Nick Henscheid
% Date:     6-2019
% Info:     Class definition file for the abstract SupportSet object.  All
%           support sets must inherit from SupportSet and are thus required
%           to implement the methods outlined in this class.
% Inputs:   
% Notes: 
% To Do:
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

classdef (Abstract) SupportSet < handle
    properties (Access=private)
        dim;
    end
    
    methods (Abstract)
        plot(obj);           % Must define geometry plotting function
        BdyFun(obj,theta);   % Define boundary function (for spherical coord)
        IsBdy(obj,r);        % Boundary checking function
        IsBdyPlus(obj,r,s);  % Outflow bdy checking
        IsBdyMinus(obj,r,s); % Inflow bdy checking
        ExitTimes(obj,r,s);  % Exit time calculation 
        UnifGrid(obj,N);     % Returns a uniform grid for the support set
        UnifSample(obj,N);   % Returns N uniformly random samples from the support set
    end
end