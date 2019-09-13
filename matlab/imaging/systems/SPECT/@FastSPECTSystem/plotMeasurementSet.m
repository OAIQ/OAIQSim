 function plotMeasurementSet(obj,gamma,j)
            % This function plots the portion of the outflow boundary that
            % is accessable given the pinhole configuration 
            figure;
            obj.plotSystem;
            Beta = zeros(2,length(j));
            bfun = obj.bdyfun;
            for i=1:length(j)
                beta = obj.dets(j(i)).getVisibleLimits(obj.bdyfun);
                Beta(1,i) = min(beta(:));
                Beta(2,i) = max(beta(:));
            end
            for l=1:length(j)
                betamin = Beta(1,l);
                betamax = Beta(2,l);
                pind = obj.dets(j(l)).detd-obj.dets(j(l)).pind;
                npin = obj.dets(j(l)).npin;
                detalpha = obj.dets(j(l)).alpha;
                pinw = obj.dets(j(l)).pinw;
                pinc = obj.dets(j(l)).pinc;
                theta = [cos(detalpha),sin(detalpha)];
                thetap = [-sin(detalpha),cos(detalpha)];
                
                for i=1:length(gamma)
                    if(obj.isAngleBetween(betamin,betamax,gamma(i)))
                        x  = [bfun(gamma(i))*cos(gamma(i)),bfun(gamma(i))*sin(gamma(i))];
                        for k=1:npin
                            pl = pind*theta + (pinc(k)-pinw/2)*thetap;
                            pr = pind*theta + (pinc(k)+pinw/2)*thetap;
                            X = [x;pl;pr];
                            fill(X(:,1),X(:,2),'b');alpha(0.2);
                        end
                    end
                end
            end
            
            
 end