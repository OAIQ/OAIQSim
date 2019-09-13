function plot(obj)
    if(obj.dim == 2)
        l = obj.detl(1); d = obj.detpos(1); alph = obj.detpos(2);
        w = obj.pinw; p = obj.pind;
        theta = [cos(alph),sin(alph)]; thetap = [-sin(alph),cos(alph)];
        xd1 = d*theta(1) - (l/2)*thetap(1); yd1 = d*theta(2) - (l/2)*thetap(2);
        xd2 = d*theta(1) + (l/2)*thetap(1); yd2 = d*theta(2) + (l/2)*thetap(2);
        xp  = zeros(obj.npin+1,2); yp  = zeros(obj.npin+1,2);
        xp(1,1) = (d-p)*theta(1) - (l/2)*thetap(1); yp(1,1) = (d-p)*theta(2) - (l/2)*thetap(2);
        xp(end,2) = (d-p)*theta(1) + (l/2)*thetap(1); yp(end,2) = (d-p)*theta(2) + (l/2)*thetap(2);
        for i=1:obj.npin
            xp(i,2)   = (d-p)*theta(1) + (obj.pinc(i)-w/2)*thetap(1);
            xp(i+1,1)   = (d-p)*theta(1) + (obj.pinc(i)+w/2)*thetap(1);
            yp(i,2)   = (d-p)*theta(2) + (obj.pinc(i)-w/2)*thetap(2);
            yp(i+1,1)   = (d-p)*theta(2) + (obj.pinc(i)+w/2)*thetap(2);
        end

        line([[xd1;xd2],xp'],[[yd1;yd2],yp'],'LineWidth',2,'Color','black');%xlim([-2.5,2.5]);ylim([-2.5,2.5]);
    end
    if(obj.dim ==3)
        % Center point is at coordinate [ax ay az].
        if(strcmp(obj.coordsys,'cart'))
            ax = obj.detpos(1);  ay = obj.detpos(2);  az = obj.detpos(3);
        elseif(strcmp(obj.coordsys,'cyl'))
            % First draw scintillator 
            drawrectcyl(obj.detpos,obj.detl);
            % Then draw pinhole
            pinholepos = obj.detpos;
            pinholepos(1) = pinholepos(1) - obj.pind;
            hold on;
            drawrectcyl(pinholepos,obj.detl);     
        elseif(strcmp(obj.coordsys,'sph'))
            [ax,ay,az] = sph2cart(obj.detpos(2),obj.detpos(3),obj.detpos(1));
        end
        % Full-width of each side of cube.
    end
end

function drawrectcyl(coords,l)
    % This function creates a planar rectangle in 3D, positioned in cylindrical
    % coordinates and with size l = [lw,lh,ld]

    rho = coords(1);
    phi = coords(2);
    z   = coords(3);

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
         ];
    patch_args = { 'FaceColor', 'b', 'FaceAlpha', 0.3 };  

    patch( 'XData', X(1,:), 'YData', X(3,:), 'ZData', X(2,:), patch_args{:} )
    %Note: Y and Z are exchanged for patient axis coords

end



