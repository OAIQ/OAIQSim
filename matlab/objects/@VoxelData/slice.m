function slice(obj)
            if(obj.dim>2)
                fig = figure('Position',[500,500,1600,1600]);
                panel1 = uipanel('Units','pixels','Position',[10,65,1600,1320]);
                fprintf('Rendering isosurface...')
                hpatch = patch(isosurface(obj.X(:,:,:,1),obj.X(:,:,:,2),obj.X(:,:,:,3),obj.Y,(min(obj.Y(:))+max(obj.Y(:)))/2));
                isonormals(obj.X(:,:,:,1),obj.X(:,:,:,2),obj.X(:,:,:,3),obj.Y,hpatch)
                hpatch.FaceColor = 'k';
                hpatch.FaceAlpha = 0.15;
                hpatch.EdgeColor = 'none';
                fprintf('done!\n');
                hold on;
                h = slice(obj.X(:,:,:,1),obj.X(:,:,:,2),obj.X(:,:,:,3),obj.Y,[],[],0.5,'cubic');
                hold off;
                set(gca,'Parent',panel1);
                set(h,'EdgeColor','none');
                axis image;
                xlim([0,1]);ylim([0,1]);zlim([0,1]);
                slx = uicontrol('Style','slider','Position',[10,35,100,30],'Value',0.5);
                xtxt = uicontrol('Style','text','String',['x= ',num2str(slx.Value)],'Position',[10,10,100,30]);
                sly = uicontrol('Style','slider','Position',[120,35,100,30],'Value',0.5);
                ytxt = uicontrol('Style','text','String',['y= ',num2str(sly.Value)],'Position',[120,10,100,30]);
                slz = uicontrol('Style','slider','Position',[230,35,100,20],'Value',0.5);
                ztxt = uicontrol('Style','text','String',['z= ',num2str(slz.Value)],'Position',[230,10,100,20]);
                slnx = uicontrol('Style','slider','Position',[340,35,100,20],'Min',-1,'Max',1);
                nxtxt = uicontrol('Style','text','String',['nx= ',num2str(slnx.Value)],'Position',[340,10,100,20]);      
                slny = uicontrol('Style','slider','Position',[450,35,100,20],'Min',-1,'Max',1);
                nytxt = uicontrol('Style','text','String',['ny= ',num2str(slny.Value)],'Position',[450,10,100,20]);
                slnz = uicontrol('Style','slider','Position',[560,35,100,20],'Value',1,'Min',-1,'Max',1);
                nztxt = uicontrol('Style','text','String',['nz= ',num2str(slnz.Value)],'Position',[560,10,100,20]);
                isolev = uicontrol('Style','slider','Position',[670,35,100,20],'Value',(min(obj.Y(:))+max(obj.Y(:)))/2,'Min',min(obj.Y(:)),'Max',max(obj.Y(:)));
                isotxt = uicontrol('Style','text','String',['isolevel= ',num2str(isolev.Value)],'Position',[670,10,150,20]);
                
                clims = caxis;
                nplottxt = uicontrol('Style','text','String','nplot = ','Position',[10,735,50,20]);
                nplot    = uicontrol('Style','edit','String',256,'Position',[70,735,100,20]);
                cmintxt  = uicontrol('Style','text','String','cmin = ','Position',[180,735,50,20]);
                cmin     = uicontrol('Style','edit','String',num2str(clims(1)),'Position',[240,735,100,20]);
                cmaxtxt  = uicontrol('Style','text','String','cmax = ','Position',[350,735,50,20]);
                cmax     = uicontrol('Style','edit','String',num2str(clims(2)),'Position',[410,735,100,20]);
                
                handles = struct('slx',slx,'sly',sly,'slz',slz,'slnx',slnx,'slny',slny,'slnz',slnz,'isolev',isolev,...
                                'xtxt',xtxt,'ytxt',ytxt,'ztxt',ztxt,'nxtxt',nxtxt,'nytxt',nytxt,'nztxt',nztxt,...
                                'isotxt',isotxt,'nplot',nplot,'cmin',cmin,'cmax',cmax);
                set(slx,'Callback',{@obj.slicecallback,panel1,handles})
                set(sly,'Callback',{@obj.slicecallback,panel1,handles})
                set(slz,'Callback',{@obj.slicecallback,panel1,handles})
                set(slnx,'Callback',{@obj.slicecallback,panel1,handles})
                set(slny,'Callback',{@obj.slicecallback,panel1,handles})
                set(slnz,'Callback',{@obj.slicecallback,panel1,handles})
                set(nplot,'Callback',{@obj.slicecallback,panel1,handles})
                set(cmin,'Callback',{@obj.slicecallback,panel1,handles})
                set(cmax,'Callback',{@obj.slicecallback,panel1,handles})
                set(isolev,'Callback',{@obj.isocallback,panel1,handles})
                
                hold on; 
                quiver3(slx.Value,sly.Value,slz.Value,0,0,0.2,'LineWidth',1.5,'Marker','.','MarkerSize',12,'Color','k');
                hold off;
            end
            % plot image
            xplot = obj.X(:,1,1);
            if(obj.dim==2)
                if(isreal(obj.Y))
                    imagesc(xplot,xplot,obj.Y');axis image;colorbar;
                    xlabel('$x$ (cm)','Interpreter','latex');
                    ylabel('$y$ (cm)','Interpreter','latex');
                    set(gca,'FontSize',18);
                else
                    figure('Position',[100,100,1600,800]);
                    subplot(1,2,1);imagesc(xplot,xplot,real(obj.Y'));axis image;title('Real Part');colorbar;
                    xlabel('$x$ (cm)','Interpreter','latex');
                    ylabel('$y$ (cm)','Interpreter','latex');
                    set(gca,'FontSize',18);
                    subplot(1,2,2);imagesc(xplot,xplot,imag(obj.Y'));axis image;title('Imaginary Part');colorbar;
                    xlabel('$x$ (cm)','Interpreter','latex');
                    ylabel('$y$ (cm)','Interpreter','latex');
                    set(gca,'FontSize',18);
                end
                % set y-axis to correct orientation
                set(gca, 'YDir', 'normal');
            end
end



