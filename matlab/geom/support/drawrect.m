function drawrect(rho,phi,z,l)
% This function creates a planar rectangle in 3D, positioned in cylindrical
% coordinates and with size l = [lw,lh,ld]

theta  = [cos(phi);sin(phi);0];
thetap = [-sin(phi);cos(phi);0];
khat   = [0;0;1];

x0 = rho*theta + [0;0;z];

lw = l(1);
lh = l(2);
ld = l(3);

hw = lw/2;
hh = lh/2;
hd = ld/2;


X = [x0 + hw*khat - hh*thetap,...
     x0 - hw*khat - hh*thetap,...
     x0 - hw*khat + hh*thetap,...
     x0 + hw*khat + hh*thetap
     ]
patch_args = { 'FaceColor', 'b', 'FaceAlpha', 0.3 };  

patch( 'XData', X(1,:), 'YData', X(3,:), 'ZData', X(2,:), patch_args{:} )

% Flip Y and Z for patient axis coords

end


