function p = ParseInputs(~,varargin)
    p = inputParser;
    p.CaseSensitive = true;
    
    det_type_default = 'pp';
    dim_default = 2;
    num_det_default = 1;
    issupp = @(f)isa(f,'SupportSet');
    
    addParameter(p,'dim',dim_default);
    addParameter(p,'num_det',num_det_default);
    parse(p,varargin{:});
    if(p.Results.dim==2)
        num_det = p.Results.num_det;
        dt = 2*pi/num_det;
        th = 0:dt:(2*pi-dt);
        for i=1:num_det
            dets_default(i)  = PPDetector('dim',2,'detpos',[12,th(i)],'detl',[4,1]);
        end
        addParameter(p,'dets',  dets_default);
        support_default  = CircleSupport(10); % 10cm radius
        addParameter(p,'support',support_default,issupp);
    elseif(p.Results.dim==3)
        num_det = p.Results.num_det;
        dt = 2*pi/num_det;
        th = 0:dt:(2*pi-dt);
        for i=1:num_det
            dets_default(i)  = PPDetector('dim',3,'detpos',[12,th(i),5],'detl',[4,4,1]);
        end
        addParameter(p,'dets',  dets_default);
        support_default  = CylinderSupport(10,10);
        addParameter(p,'support',support_default,issupp);
    end
    addParameter(p,'det_type', det_type_default);
    
    parse(p,varargin{:});
end