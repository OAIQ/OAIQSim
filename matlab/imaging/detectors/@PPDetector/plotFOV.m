function plotFOV(obj,x,support)
    plot(obj); hold on;
    alp = obj.detpos(2);
    for i=1:length(x)
        beta = obj.getPinholeAngles(x(i));
        theta  = [cos(alp);sin(alp)];
        thetap = [-sin(alp);cos(alp)];
        p = obj.detpos(1)*theta + x(i)*thetap;
        plot(p(1),p(2),'.','MarkerSize',10);
        for j=1:obj.npin
            x1 = @(t) p(1) - t*cos(beta(1,j));
            y1 = @(t) p(2) - t*sin(beta(1,j));
            x2 = @(t) p(1) - t*cos(beta(2,j));
            y2 = @(t) p(2) - t*sin(beta(2,j));
            tau1 = support.ExitTimes([p(1),p(2)],-[cos(beta(1,j)),sin(beta(1,j))]);
            tau2 = support.ExitTimes([p(1),p(2)],-[cos(beta(2,j)),sin(beta(2,j))]);
            tplot = linspace(min(tau1(1),tau1(2)),0,2);
            X1 = x1(tplot);Y1 = y1(tplot);
            tplot = linspace(min(tau2(1),tau2(2)),0,2);
            X2 = x2(tplot);Y2 = y2(tplot);
            plot(X1,Y1);
            plot(X2,Y2);
            X = unique([[X1(:),Y1(:)];[X2(:),Y2(:)]],'rows');

            fill(X(:,1),X(:,2),'b','LineStyle','none'); alpha(0.2);
        end
    end
end