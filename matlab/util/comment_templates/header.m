% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
% File:     RTELumpy.m
% Author:   Nick Henscheid
% Date:     6-2019
% Info:     Class definition file for the RTELumpy object, a ballistic RTE
%           solver for lumpy background objects.
% Inputs:   An object of type RTELumpy has the following optional inputs:
%           - 'Xi' a source function of LumpyBgnd type
%           - 'Mu' an attenuation function of LumpyBgnd type
%           - 'cm' a medium speed of light
%           - 'hmin' a minimum step length for the ballistic RTE solver 
%           - Defaults will be chosen for each if not specified by user.
% Notes:    - The dimension is fixed by the Xi and Mu functions (which must
%             obviously agree)
% To Do:    - 
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------