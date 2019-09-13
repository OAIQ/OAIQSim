% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:    test_lumpy.m
% Author:  Nick Henscheid
% Date:    3-2017
% Info:    test_lumpy is a unit test for the LumpyBgnd class.
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

function TestVoxel(test_settings)


% Test basic object creation

GPU = test_settings.GPU;
DISPLAY = test_settings.DISPLAY;


V = VoxelData;

if(DISPLAY)
    plot(V);
end







end % TestVoxel



 
