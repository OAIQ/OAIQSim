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

function TestLumpy(test_settings)

nplotx = 4;
nploty = 4;

% Test basic object creation

GPU = test_settings.GPU;
DISPLAY = test_settings.DISPLAY;


disp('Creating default lumpy background object')
L = LumpyBgnd  % Default object: 2D lumpy background with default parameters and randomly generated centers

L.gpu = GPU;

if(DISPLAY)
    fig_window = figure;
    fig_window.Position = [fig_window.Position(1),fig_window.Position(2),800,800];
end

% Test plot function 
if(DISPLAY)
    fprintf('Testing plot function\n');
    ax = subplot(nplotx,nploty,1);
    plot(L,ax);title('Result of L.plot()');colorbar;
end

% Test 2D evaluation

fprintf('Testing basic 2D evaluation on the default grid (%dx%d), returning simple array\n',L.N,L.N) 

tic
u = L.eval();
t = toc;

fprintf('Success! Time Elapsed = %fs\n',t);

if(DISPLAY)
    subplot(nplotx,nploty,2);
    imagesc(u);axis image;set(gca,'YDir','normal');title('Result of L.eval()');colorbar;
end


% Test randomization 

fprintf('Randomizing lump centers\n');
tic
L.randomize();
fprintf('Success! Time elapsed = %fs\n',toc);

if(DISPLAY)
    ax = subplot(nplotx,nploty,3);
    plot(L,ax);title('Re-randomized'); colorbar;
end


% Test grid refinement 

fprintf('Refining the evaluation grid to %dx%d & re-evaluating\n',512,512);
L.N = 512;
tic
u   = eval(L);
fprintf('Success! Time elapsed = %f\n',toc);

if(DISPLAY)
    subplot(nplotx,nploty,4);
    imagesc(u);axis image;set(gca,'YDir','normal');title('Result of L.eval() w/ L.N = 512');colorbar;
end


% Test lump size change
fprintf('Changing the lump size to 0.001 (isotropic) and re-evaluating.  Note the warning!\n');
L.cov = 0.001;
tic
u   = eval(L);
fprintf('Success! Time elapsed = %f\n',toc);


if(DISPLAY)
    ax = subplot(nplotx,nploty,5);
    plot(L,ax);title(sprintf('Lump radius = %1.1e',L.cov(1)));colorbar;
end

% Turn off warnings and change the lump size again
fprintf('Turning off warnings with L.TurnOffWarnings()\n');
L.TurnOffWarnings;
L.cov = 0.01;
tic
u   = eval(L);
fprintf('Success! Time elapsed = %f\n',toc);


if(DISPLAY)
    ax = subplot(nplotx,nploty,6);
    plot(L,ax);title(sprintf('Lumpy radius = %1.1e',L.cov(1)));colorbar;
end


% Test anisotropic-but-diagonal lump
fprintf('Changing the lump covariance matrix to anisotropic\n');
L.cov = diag([0.01,0.001]);
tic
u   = eval(L);
fprintf('Success! Time elapsed = %f\n',toc);

if(DISPLAY)
    ax = subplot(nplotx,nploty,7);
    plot(L,ax);title('Anisotropic Covariance');colorbar;
end


% Test rotated anisotropic lump
fprintf('Rotating the lump covariance matrix\n');
theta = pi/4;
R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
L.cov = R*L.cov*R';
tic
u   = eval(L);
fprintf('Success! Time elapsed = %f\n',toc);

if(DISPLAY)
    ax = subplot(nplotx,nploty,8);
    plot(L,ax);title('Rotated Covariance by $\pi/4$');colorbar;
end

% Test non-uniform anisotropic lumps
fprintf('Testing non-uniform lumps (random rotation)\n');
for i=1:L.K-1
    theta = 2*pi*rand();
    R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    L.cov = [L.cov,R*L.cov(:,2*i-1:2*i)*R'];
end

tic
u   = eval(L);
fprintf('Success! Time elapsed = %f\n',toc);

if(DISPLAY)
    ax = subplot(nplotx,nploty,9);
    plot(L,ax);title('Lumps randomly rotated');colorbar;
end


% Test non-uniform lump amplitude 
fprintf('Testing non-uniform lump amplitude (lognormal distribution)\n');

L.b0 = L.b0*exp(randn(L.K,1));

tic
u   = eval(L);
fprintf('Success! Time elapsed = %f\n',toc);

if(DISPLAY)
    ax = subplot(nplotx,nploty,10);
    plot(L,ax);title('Lumps randomly scaled (lognormal distr.)');colorbar;
end

% Test changing mean number of lumps




% Test 3D evaluation 



 
