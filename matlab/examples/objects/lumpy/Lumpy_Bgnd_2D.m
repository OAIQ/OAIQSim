% -------------------------------------------------------------------------
% ********************OAIQSIM SIMULATION TOOLBOX***************************
%
% File:     Lumpy_Bgnd_2D.m
% Author:   Nick Henscheid
% Date:     4-2019,9-2019
% Info:     This script demonstrates the lumpy background object functionality in 2D
% Notes:   
% To Do: 
%
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

close all;
% NOTE: This script will produce a figure similar to figure XXX in the User's Guide, but the  
% realizations will be difference since LumpyBgnd.randomize calls
% rng('shuffle') every time.
%
% First create the lumpy background object with a rectangular support set 
L = LumpyBgnd('support',RectSupport('L',[0,1;0,1]));  
% Adjust some settings
L.TurnOffWarnings;  % Don't show warning messages (comment out if you want to see them).
L.Kbar = 200;       % Adjust the mean number of lumps
L.N = 256;          % Adjust the mesh size for plotting
L.gpu = 0;          % Run on the CPU 

% Display a realization 
plot(L);title('Realization of LumpyBgnd');

%% Randomize the field and re-plot
L.randomize;
fig2 = figure;
subplot(2,4,1);
plot(L);title('L.cov = 0.005','FontSize',12);

% Change the lump width, but keep it isotropic and keep the same centers
% Note: if line 25 is commented out, a warning will display in the command window, informing  
% you that the original lump centers were selected based on the original
% lump width and hence boundary artifacts might occur.

L.cov = 0.001;
subplot(2,4,2);
plot(L);title('L.cov = 1e-3','FontSize',12);

L.cov = 0.0001;
subplot(2,4,3);
plot(L);title('L.cov = 1e-4','FontSize',12);

L.cov = 0.1;
subplot(2,4,4);
plot(L);title('L.cov = 0.1','FontSize',12);

%% Change the lump to be anisotropic (first vertical then horizontal)
L.cov = diag([0.001,0.01]);
subplot(2,4,5);
plot(L);title('L.cov = diag([1e-3,1e-2])','FontSize',12);

L.cov = diag([0.05,0.0005]);
subplot(2,4,6);
plot(L);title('L.cov = diag([5e-2,5e-4])','FontSize',12);
%% Rotate lumps by a fixed rotation matrix 
theta = pi/4;  % Rotation angle
var1  = 0.01;   % var (rotated x-hat)
var2  = 0.001; % var (rotated y-hat)

R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
K = R*diag([var1,var2])*R';

L.cov = K;

subplot(2,4,7);
plot(L);title(sprintf('L.cov = [%1.1e,%1.1e;%1.1e,%1.1e]',L.cov(1,1),L.cov(1,2),L.cov(2,1),L.cov(2,2)),'FontSize',12);
%% Make the lumps non-uniform in size (randomly generated isotropic covariance)
lump_mean = log(0.07);
lump_std  = log(0.12);
L.cov = exp(lump_mean+lump_std*randn())*eye(2);  % Lognormally distributed
L.b0  = 1/sqrt((2*pi)^2*det(L.cov));  %PDF/Gaussian mixture model scaling
for i=1:L.K-1
    covtemp = exp(lump_mean+lump_std*randn())*eye(2);
    L.cov = [L.cov,covtemp];
    L.b0  = [L.b0,1/sqrt((2*pi)^2*det(covtemp))];
end

subplot(2,4,8);
plot(L);title('Indepdent cov for ea. lump','FontSize',12);


