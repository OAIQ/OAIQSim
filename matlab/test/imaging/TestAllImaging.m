% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:     TestAllImaging.m
% Author:   Nick Henscheid
% Date:     6-2019
% Info:     Script to run all imaging tests
%               
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

disp('Testing All Imaging Methods');

GPU = sign(gpuDeviceCount);
DISPLAY = usejava('desktop');

test_settings.GPU = GPU;
test_settings.DISPLAY = DISPLAY;

TestFastSPECT(test_settings);
