% Test RTELumpy 
close all

% Support is a unit cube 
S  = RectSupport('L',[0,1;0,1;0,1]);
% Define the source and attenuation functions
Xi = LumpyBgnd('support',S,'Kbar',30,'cov',0.001);
Mu = LumpyBgnd('support',S,'Kbar',20,'cov',0.1);
% Define the RTE Solver
W = RTELumpy('Xi',Xi,'Mu',Mu);

% Set the step length for the RTE solver
W.hmin = 1e-3;

figure; subplot(1,2,1); plot(Xi,gca); title('Isosurface of source');
subplot(1,2,2); plot(Mu,gca); title('Isosurface of attenuation');
% Assume zero attenuation at first
Mu.b0 = 0.0;

%%
% Assume the detector is placed "above" the support and we measure
% collimated rays pointing in the positive z direction
nx = 128;
xtemp = linspace(0,1,nx);
[xx,yy] = meshgrid(xtemp);
R = [xx(:)';yy(:)';ones(1,nx^2)];
S = repmat([0;0;1],[1,nx^2]);

subplot(1,2,1);hold on;
plot3(R(1,:),R(2,:),R(3,:),'.');
quiver3(R(1,:),R(2,:),R(3,:),S(1,:),S(2,:),S(3,:));
hold off;
R = R(:);
S = S(:);

fprintf('Computing solution without attenuation\n');
y0 = W.eval(R,S,0);

y0 = real(reshape(y0,[nx,nx]));
fig1 = figure; subplot(1,3,1);
imagesc(xtemp,xtemp,y0,[0,max(y0(:))]);
set(gca,'YDir','normal');colorbar;
title('Solution without attenuation');
%%
fprintf('Computing voxelized solution for comparison\n');
Xi.N = 128;
dx = 1/(Xi.N-1);
u = Xi.eval;
us = sum(u.data_array,3)*dx/W.cm;
subplot(1,3,2);
xtemp2 = linspace(0,1,Xi.N);
imagesc(xtemp2,xtemp2,us,[0,max(y0(:))]);
set(gca,'YDir','normal');colorbar;
title('Voxelized solution');
%%
idx = 8:8:32;
figure;
p1=plot(xtemp,y0(:,idx));
hold on; 
p2 = plot(xtemp2,us(:,idx),'o');
title(sprintf('Comparison of RTE solution to voxelized solution ($h_{RTE}$ = %1.1e, $h_{vox} = $%1.1e)',W.hmin,1/Xi.N));
legend([p1(1),p2(1)],'RTE Solution','Voxelized solution');

%%
% Now, with attenuation.

fprintf('Computing solution with attenuation\n');
W.Mu.b0 = 1;
y1 = real(reshape(W.eval(R,S,0),[nx,nx]));
figure(fig1);
subplot(1,3,3);
imagesc(xtemp,xtemp,y1,[0,max(y0(:))]);
set(gca,'YDir','normal');colorbar;
title('Solution with attenuation');


