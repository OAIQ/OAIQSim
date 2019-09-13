% ---------------------------------------------------------------------------
% **********************CGRI SIMULATION TOOLBOX***************************
% test_lumpy is a testing script for the LumpyBgnd class.

% File: test_lumpy.m
% Author:  Nick Henscheid
% Date:    3-2017
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without technical
% support, and with no warranty, express or implied, as to its usefulness for
% any purpose.
% ---------------------------------------------------------------------------


% Test basic object creation and properties

disp('Creating default lumpy background object')
u0 = LumpyBgnd  % Default object: 2D lumpy background with default parameters and randomly generated centers


a0 = u0.eval;
imagesc(a0);title('Texture evaluated on the default grid')

nx = floor(input('Enter a grid size nx to re-evaluate the texture: '));
[x,y] = meshgrid(linspace(0,1,nx));
X = [x(:) y(:)];
a1 = u0.eval(X,[nx,nx]);
subplot(2,2,2);
imagesc(a1);title('Texture evaluated on a user-supplied grid');

s=1;
while(s==1)
    s = input('Re-randomize centers? (1 = yes, 0 = no) ');
    if(s==1)
        u0.randomize;
        a2 = u0.eval(X,[nx,nx]);
        subplot(2,2,3);
        imagesc(a2);title('Texture with re-randomized centers evaluated on the same grid')
    end
end

s=1;
while(s==1)
    r = input('Enter a lump width or 0 to exit: ');
    if(r>0)
        u0.cov = r;
        a3 = u0.eval(X,[nx,nx]);
        subplot(2,2,4);
        imagesc(a3);title('Same texture with smaller lump width')
    else
        s=0;
    end
end

