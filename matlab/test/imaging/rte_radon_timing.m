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

% Plotting?
plots = 1;


% Set up geometry
alpha = 0;
s0  = [cos(alpha),sin(alpha)];
s0p = [-sin(alpha),cos(alpha)];
beta  = 0;

nps = 10.^(2:5);
nsteps = floor(linspace(1e1,1e4,20));
hmins = linspace(1e-4,0.1,20);

nnps = length(nps);
nnsteps = length(nsteps);

times = zeros(nnps,nnsteps);
errs2 = zeros(nnps,nnsteps);
errsinf = zeros(nnps,nnsteps);

p = linspace(-1,1,nps(1));
R = [sqrt(1-p.^2);p];
S = repmat([1;0],[1,nps(1)]);
E = zeros(length(R(:)),1);

gpu = 1;
k = [5,5];

%gpu = input('Use GPU? (enter 0 or 1) ');
T = RTEFourier('k',k,'gpu',logical(gpu));


if(~usejava('desktop'))
    disp('Running in text mode...');
    plots = 0;
end

c = 299792458;
% Exact Radon transform of a Fourier mode.
phi = @(p) 2*sqrt(1-p.^2).*exp(2*pi*1i*p*dot(k,s0p)).*sinc(2*sqrt(1-p.^2)*dot(k,s0))/(4*pi*c);

% Warm-up
disp('Warming up!')
w = T.eval(R(:),S(:),E);
disp('Done warming up.');


for i=1:nnps
    for j=1:nnsteps
        progressbar(i/nnps,j/nnsteps)
        nstep = nsteps(j);
        T.hmin = hmins(nnsteps-(j-1));
        np = nps(i);
        p = linspace(-1,1,np);
        R = [sqrt(1-p.^2);p];
        S = repmat([1;0],[1,np]);

        Rin = R(:);
        Sin = S(:);
        E = zeros(length(Rin),1);
        tic
        w = T.eval(Rin,Sin,E);
        times(i,j) = toc;
        w_e = phi(p).';
        errs2(i,j) = norm(w-w_e)/norm(w_e);
        errsinf(i,j) = norm(w-w_e,inf)/norm(w_e,inf);
        fprintf('Evaluation took %f seconds. (err2, errinf) = (%2.6e,%2.6e)\n',times(i,j),errs2(i,j),errsinf(i,j));
    end
end
%
%figure;


if(plots)
    subplot(2,2,4);
    plot(p,real(w),'linewidth',1.5);hold on;
    plot(p,imag(w),'linewidth',1.5);hold on;
    plot(p,real(w_e),'--','linewidth',1.5)
    plot(p,imag(w_e),'--','linewidth',1.5)
    title(sprintf('\\makebox[4.5in][c]{Numerical and exact 2D Radon transform of $\\exp(2\\pi i(%1.0fx+%1.0fy))/4\\pi$}\n\\makebox[4.5in][c]{error(re,im) = (%2.3e,%2.3e)}',k(1),k(2),norm(real(w-w_e))/norm(real(w_e)),norm(imag(w-w_e))/norm(imag(w_e))),'interpreter','latex','fontsize',14);
    leg = legend('numerical (re)','numerical (im)','exact (re)','exact (im)');set(leg,'fontsize',14);
else
    disp('Times:')
    disp(times)
    disp('Errs:')
    disp(errs2)
    disp(errsinf)
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


