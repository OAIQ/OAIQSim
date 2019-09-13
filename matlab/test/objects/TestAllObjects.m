% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:     TestAllObjects.m
% Author:   Nick Henscheid
% Date:     6-2019
% Info:     Script to run all object tests
%                       
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

close all;


disp('Testing All Object Methods');


GPU = sign(gpuDeviceCount);
DISPLAY = usejava('desktop');

test_settings.GPU = GPU;
test_settings.DISPLAY = DISPLAY;


% Lumpy Background test
%TestLumpy(test_settings);



% VoxelData Test
TestVoxel(test_settings);
