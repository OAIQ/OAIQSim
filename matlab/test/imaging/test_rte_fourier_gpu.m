% A test for the RTEFourier object

close all


%% 3D, evaluation points lie on surface of sphere.

% Number of polar and azimuth angles
nth = 256;
nph = 256;

dth = pi/(nth+1);
th  = dth:dth:(pi-dth);
dph = pi/(nph+1);
ph  = (-pi/2+dph):dph:(pi/2-dph);
[thg,phg] = ndgrid(th,ph);
thg = thg(:)';
phg = phg(:)';
x = sin(thg).*cos(phg); 
y = sin(thg).*sin(phg);
z = cos(thg);
xx = reshape(x,[nth,nph]);
yy = reshape(y,[nth,nph]);
zz = reshape(z,[nth,nph]);
R = [x;y;z];
% Directions all point to the "right" (along xhat direction)
S = repmat([1;0;0],[1,nth*nph]);

%%
figure; plot3(x(:),y(:),z(:),'.');axis equal;
xlabel('$x$','Interpreter','latex');ylabel('$y$','Interpreter','latex');zlabel('$z$','Interpreter','latex');
hold on; quiver3(x(:),y(:),z(:),S(1,:)',S(2,:)',S(3,:)'); set(gca,'Clipping','off');
title('Evaluation points and directions','FontSize',20);

%%
k = [1,0,0];
f = @(x,y,z) exp(2*pi*i*(k(1)*x+k(2)*y+k(3)*z));

fg = f(xx,yy,zz);

W = RTEFourier('k',k,'cm',1,'gpu',true);
w = W.eval(R(:),S(:),0.0);
wre = reshape(real(w),[nth,nph]);
wim = reshape(imag(w),[nth,nph]);

%%
figure;
subplot(1,2,1); imagesc(th,ph,wre); axis image; set(gca,'YDir','normal'); title('Real part of $w(r,\hat{s})$'); colorbar; 
subplot(1,2,2); imagesc(th,ph,wim); axis image; set(gca,'YDir','normal'); title('Imaginary part of $w(r,\hat{s})$'); colorbar;

figure; 
subplot(2,2,1);
surf(xx,yy,zz,wre,'EdgeColor','none');axis equal; title('Real part of $w(r,\hat{s})$'); colorbar;
subplot(2,2,2);
surf(xx,yy,zz,wim,'EdgeColor','none');axis equal; title('Imaginary part of $w(r,\hat{s})$'); colorbar;

subplot(2,2,3);
surf(xx,yy,zz,real(fg),'EdgeColor','none');axis equal; title('Real part of $f(r)$'); colorbar;
subplot(2,2,4); 
surf(xx,yy,zz,imag(fg),'EdgeColor','none');axis equal; title('Imaginary part of $f(r)$'); colorbar;