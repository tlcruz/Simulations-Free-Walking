function [str] = GetStraightness(dt, locSpk, params)
actst = dt.actSt;
actst(dt.WD < params.mDistWall) = 0;
addpath('StraightnessFunctions');

locsF = locSpk;
pksIF = 7*ones(length(locSpk),1);
pksEF = 7*ones(length(locSpk),1);
[FBouts] = GetFbouts(actst, dt.Vf, locsF, pksIF,pksEF, params);

STR = [];
ANGD = [];
DST = [];
N = [];

for pp = 1 : length(FBouts)
    % Calculate Straightness
    if length(FBouts{pp}) > max(params.windStr+1, params.minStrB);
        xaux = dt.X(FBouts{pp});
        yaux = dt.Y(FBouts{pp});
        vtaux = dt.Vt(FBouts{pp});
        vraux = dt.Vr(FBouts{pp});
        nv = length(xaux) - params.windStr;
        dm = [];
        ss = [];
        for l = 1 : nv
            pi = [xaux(l),yaux(l),0];
            pf = [xaux(l+params.windStr),yaux(l+params.windStr),0];
            dis = sum(vtaux((l):(l+params.windStr)))/60;
            pt = [xaux(l+floor(params.windStr/2)),...
                yaux(l+floor(params.windStr/2)),0];
            dr = point_to_line(pt,pi,pf);
            dm = vertcat(dm, dr);
            ss = vertcat(ss, dis);
        end
        STR = vertcat(STR, sum(ss)/sum(dm));
        DST = vertcat(DST, sum(ss));
        N = vertcat(N, length(xaux));
%         ANGD = vertcat(ANGD, sum(abs(vraux))/sum(vtaux));
        ANGD = vertcat(ANGD, mean(abs(vraux)));
    end
end

str.STR = STR;
str.DST = DST;
str.ANGD = ANGD;
str.N = N;
str.locsF = locsF;
str.pksIF = pksIF;
str.pksEF = pksEF;
str.FBouts = FBouts;

end