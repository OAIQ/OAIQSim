function beta = getPinholeAngles(obj,x)
    l = obj.detl(1);
    if(abs(x)>l/2)
        error('x is out of bounds!')
    end
    w = obj.pinw;
    d = obj.pind;
    c = obj.pinc;
    alpha = obj.detpos(2);

    beta = zeros(2,obj.npin);

    for i=1:obj.npin
        if(x>c(i)+w/2)
            beta(1,i) = 3*pi/2 + alpha - atan(d/abs(x-(c(i)+w/2)));
        else
            beta(1,i) = pi/2 + alpha + atan(d/abs(x-(c(i)+w/2)));
        end
        if(x>c(i)-w/2)
            beta(2,i) = 3*pi/2 + alpha - atan(d/abs(x-(c(i)-w/2)));
        else
            beta(2,i) = pi/2 + alpha + atan(d/abs(x-(c(i)-w/2)));
        end
    end
end