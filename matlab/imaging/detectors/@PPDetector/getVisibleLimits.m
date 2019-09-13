function beta = getVisibleLimits(obj,bdyfun)
    if(obj.dim==2)
        beta = zeros(2,obj.npin);
        for i=1:obj.npin
            sl = obj.pinc(i)-obj.pinw/2; 
            sr = obj.pinc(i)+obj.pinw/2;
            d = obj.detl(end) - obj.pind;
            xl = d*cos(obj.detpos(2)) - sl*sin(obj.detpos(2));
            yl = d*sin(obj.detpos(2)) + sl*cos(obj.detpos(2));
            xr = d*cos(obj.detpos(2)) - sr*sin(obj.detpos(2));
            yr = d*sin(obj.detpos(2)) + sr*cos(obj.detpos(2));
            options = optimoptions('fsolve','Display','none');
            F = @(x)[bdyfun(x(1)).*cos(x(1)) - x(2).*sin(x(1))-xl,bdyfun(x(1)).*sin(x(1)) + x(2).*cos(x(1))-yl];
            sol = fsolve(F,[obj.detpos(2)-pi/2,0],options);
            beta(1,i) = sol(1);
            F = @(x)[bdyfun(x(1)).*cos(x(1)) - x(2).*sin(x(1))-xr,bdyfun(x(1)).*sin(x(1)) + x(2).*cos(x(1))-yr];
            sol = fsolve(F,[obj.detpos(2)+pi/2,0],options);     
            beta(2,i) = sol(1);
        end
    else
        error('Only works for dim = 2!');
    end
end