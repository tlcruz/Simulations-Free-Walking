function [PSS] = GetPSS(dt, locSpk, params)
actst = dt.actSt;
actst(dt.WD < params.mDistWall) = 0;
addpath('StraightnessFunctions');

locsF = locSpk;
pksIF = 7*ones(length(locSpk),1);
pksEF = 7*ones(length(locSpk),1);
[FBouts] = GetFbouts(actst, dt.Vf, locsF, pksIF,pksEF, params);
[Bouts] = GetFbouts(actst, dt.Vf, [], [],[], params);
npss = [];
pss = [];
vrl = [];
pvrl = [];
vrr = [];
nvrl = [];
nvrr = [];
pvrr = [];
for i = 1 : length(Bouts)
    if length(Bouts{i}) > max(params.windStr+1, params.minStrB);
        vrb = dt.Vr(Bouts{i});
        vfb = dt.Vf(Bouts{i});
        [BoutVr,~,~,BoutP] = GetProbVec(vrb,vfb, params.delta,1,-inf,inf);        
        npss = vertcat(npss,length(BoutP));
        pss = horzcat(pss,mean(BoutP));
        
        inds = vrb;
        inds(inds>0) = 1;
        inds(inds<=0) = -1;
        aux1 = inds(1:end-params.delta+1);
        aux2 = inds(params.delta:end);

        if ~isempty(find(aux1 == 1, 1))
            pvrr = horzcat(pvrr, ...
                length(find((aux1+aux2)/2 == 1))/(length(find(aux1 == 1))+1));
            nvrr = vertcat(nvrr, length(find(aux1 == 1)));
            vrr = horzcat(vrr, abs(mean(vrb(aux1 == 1))));
        end
        if ~isempty(find(aux1 == -1, 1))
            pvrl = horzcat(pvrl, ...
                length(find((aux1+aux2)/2 == -1))/(length(find(aux1 == -1))+1));
            nvrl = vertcat(nvrl, length(find(aux1 == -1)));
            vrl = horzcat(vrl, abs(mean(vrb(aux1 == -1))));
        end
        
    end
end
% for i = 1 : length(FBouts)
%     if length(FBouts{i}) > max(params.windStr+1, params.minStrB);
%         vrb = dt.Vr(FBouts{i});
%         vfb = dt.Vf(FBouts{i});
%         [BoutVr,~,~,BoutP] = GetProbVec(vrb,vfb, params.delta,1,-inf,inf);
%         npss = vertcat(npss,length(BoutP));
%         pss = horzcat(pss,mean(BoutP));
%         vra1 = vrb((params.delta+1):(end-20));
%         vra2 = vrb(1:(end-params.delta-20));
%         
%         inds = vrb;
%         inds(inds>0) = 1;
%         inds(inds<=0) = -1;
%         aux1 = inds(1:end-params.delta+1);
%         aux2 = inds(params.delta:end);
% 
%         if ~isempty(find(aux1 == 1, 1))
%             pvrr = horzcat(pvrr, ...
%                 length(find((aux1+aux2)/2 == 1))/(length(find(aux1 == 1))+1));
%             nvrr = vertcat(nvrr, length(find(aux1 == 1)));
%             vrr = horzcat(vrr, abs(mean(vrb(aux1 == 1))));
%         end
%         if ~isempty(find(aux1 == -1, 1))
%             pvrl = horzcat(pvrl, ...
%                 length(find((aux1+aux2)/2 == -1))/(length(find(aux1 == -1))+1));
%             nvrl = vertcat(nvrl, length(find(aux1 == -1)));
%             vrl = horzcat(vrl, abs(mean(vrb(aux1 == -1))));
%         end
%     end
% end

PSS.pss = pss*npss/sum(npss);
PSS.npss = sum(npss);
PSS.vrl = vrl*nvrl/sum(nvrl);
PSS.pvrl = pvrl*nvrl/sum(nvrl);
PSS.nvrl = sum(nvrl);
PSS.vrr = vrr*nvrr/sum(nvrr);
PSS.pvrr = pvrr*nvrr/sum(nvrr);
PSS.nvrr = sum(nvrr);
end
