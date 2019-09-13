% ---------------------------------------------------------------------------
% **********************CGRI SIMULATION TOOLBOX***************************
% File: test_rte_radon_cpu.m
% Purpose: testing the cpu rte solver with parallel beam Radon transform
%          geometry. 
% Author:  Nick Henscheid
% Date:    5-2017
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without technical
% support, and with no warranty, express or implied, as to its usefulness for
% any purpose.
% ---------------------------------------------------------------------------

% Set up geometry
alpha = 0;
s0  = [cos(alpha),sin(alpha)];
s0p = [-sin(alpha),cos(alpha)];
beta  = 0;
np = input('How many rays? ');
p = linspace(-1,1,np);
R = [sqrt(1-p.^2);p];
S = repmat([1;0],[1,np]);
k = input('Enter a wave vector in the form [k1,k2]: ');
gpu = input('Use GPU? (enter 0 or 1) ');
plots = input('Display plots? (enter 0 or 1) ');
if(~usejava('desktop'))
    disp('Running in text mode...');
    plots = 0;
end

if(plots)
    %  Plot the evaluation rays on the boundary
    tplot = linspace(0,2*pi,100);
   % figure('position',[100,100,800,800]);subplot(2,2,3);
    figure;
    plot(cos(tplot),sin(tplot));hold on;axis equal;
    plot(R(1,:),R(2,:),'.')
    quiver(R(1,:),R(2,:),S(1,:),S(2,:));
    % Plot the Fourier mode
    xplot = linspace(-1,1,256);
    [xx,yy] = meshgrid(xplot);
    suppfun = @(x,y) (x.^2+y.^2)<1;
    fmode   = @(x,y,k) exp(2*pi*1i*(k(1)*x + k(2)*y));
    %subplot(2,2,1);
    figure;
    imagesc(xplot,xplot,real(suppfun(xx,yy).*fmode(xx,yy,k))); 
    axis image; xlabel('Real part'); 
    %subplot(2,2,2);
    figure;
    imagesc(xplot,xplot,imag(suppfun(xx,yy).*fmode(xx,yy,k))); 
    axis image; xlabel('Imaginary part');
end

T = RTEFourier('k',k,'hmin',1e-4,'gpu',logical(gpu),'cm',1);

Rin = R(:);
Sin = S(:);
E = zeros(length(Rin),1);
tic
w = T.eval(Rin,Sin,E);
disp(sprintf('Evaluation took %f seconds\n',toc))
%
%figure;
%c = 299792458;
c = 1;
% Exact Radon transform of a Fourier mode.
phi = @(p) 2*sqrt(1-p.^2).*exp(2*pi*1i*p*dot(k,s0p)).*sinc(2*sqrt(1-p.^2)*dot(k,s0))/(4*pi*c);
w_e = phi(p).';
if(plots)
    %subplot(2,2,4);
    figure;plot(p,real(w),'linewidth',1.5);hold on;
    plot(p,imag(w),'linewidth',1.5);hold on;
    plot(p,real(w_e),'--','linewidth',1.5)
    plot(p,imag(w_e),'--','linewidth',1.5)
    title(sprintf('\\makebox[4.5in][c]{Numerical and exact 2D Radon transform of $\\exp(2\\pi i(%1.0fx+%1.0fy))/4\\pi$}\n\\makebox[4.5in][c]{error(re,im) = (%2.3e,%2.3e)}',k(1),k(2),norm(real(w-w_e))/norm(real(w_e)),norm(imag(w-w_e))/norm(imag(w_e))),'interpreter','latex','fontsize',14);
    leg = legend('numerical (re)','numerical (im)','exact (re)','exact (im)');set(leg,'fontsize',14);
else
    fprintf('Error (re) = %2.3e\n',norm(real(w-w_e))/norm(real(w_e)));
    fprintf('Error (im) = %2.3e\n',norm(imag(w-w_e))/norm(imag(w_e)));
end

%%  This was some Photon processing stuff - not needed here. 
% sig = 1e-2;
% %convkern = @(p) exp(-p.^2./(2*sig^2))/sqrt(2*pi*sig^2);
% h = p(2) - p(1);
% v = gauss(sig,p);
% tic
% phat = h*conv(w,v,'same');
% timeconv = toc;
% figure;
% plot(p,w);hold on;plot(p,phat);title(num2str(timeconv));


