function [outp] = StartSimulation(params)
XX = params.x0;
YY = params.y0;
TH = params.th0;
fM = 0;
fV = 0;
PER = 0;
dBS = load('DistBoutLength.mat');
locSpk = [];
nLR = 0;
N = 0;
while length(XX) < params.L
    %% Enter walking bout
    if rand > 0.99
        [rNs] = RandDist(dBS.distBoutSize, dBS.lcents, 1);
        rNs = floor(rNs);
        if params.PertOn == 1
           [per] = GetPerturbation(rNs, params);
        else
            per.on = 0;
        end
        noisevect = params.noisesize*flicker(50+rNs);
        noisevect = filter(params.noisebf,params.noiseaf,noisevect);
        %% For the duration of the walking bout
        for i = 1 : rNs
            x0 = XX(end);
            y0 = YY(end);
            ang = TH(end)+noisevect(i+20);
            vr = diff(TH);
            % Apply Control
            [ang, mF, vF] = ApplyControl(vr, ang, params, per, i);
            
            x = x0 + params.vf*cos(pi*ang/180)/60;
            y = y0 + params.vf*sin(pi*ang/180)/60;
            
            if IsSpike(XX, x, y, params)
                locSpk = vertcat(locSpk, length(XX));
                [XX, YY, TH] = SpikeTurn(x0, y0, x, y, ang, params.vf, ...
                    XX, YY, TH, params, nLR, N);
                fM = vertcat(fM(1:(end-20)), zeros(30,1));
                fV = vertcat(fV(1:(end-20)), zeros(30,1));
                PER = vertcat(PER(1:(end-20)), zeros(30,1));
            else
                nLR = nLR + sign(ang-TH(end));
                N = N + 1;
                TH = vertcat(TH, ang);
                XX = vertcat(XX, x);
                YY = vertcat(YY, y);
                fM = vertcat(fM, mF);
                fV = vertcat(fV, vF);
                if per.on == 1
                    PER = vertcat(PER, per.pert(i));
                else
                    PER = vertcat(PER, 0);
                end
            end
        end
    end
    %% Outside of walking bout
    TH = vertcat(TH, TH(end));
    XX = vertcat(XX, XX(end));
    YY = vertcat(YY, YY(end));
    fM = vertcat(fM, 0);
    fV = vertcat(fV, 0);
    PER = vertcat(PER,0);
end
[dt] = PreProcessing(XX,YY,TH);
[str] = GetStraightness(dt, locSpk, params);
[pss] = GetPSS(dt, locSpk, params);
[pertVrMat, pertMat] = ProcessPerturbations(dt, PER, params);
outp.PER = PER;
outp.pertVrMat = pertVrMat;
outp.pertMat = pertMat;
outp.strA = str;
outp.str = str.STR'*str.N/sum(str.N);
outp.angd = str.ANGD'*str.N/sum(str.N);
outp.dist = str.DST'*str.N/sum(str.N);
outp.n = sum(str.N);
outp.pss = pss.pss;
outp.npss = pss.npss;
outp.vrl = pss.vrl;
outp.pvrl = pss.pvrl;
outp.nvrl = pss.nvrl;
outp.vrr = pss.vrr;
outp.pvrr = pss.pvrr;
outp.nvrr = pss.nvrr;
outp.dt = dt;
outp.params = params;
outp.fM = fM;
outp.fV = fV;
if params.plotTraj == 1
    
    t = -0.1:0.1:2*pi;
    figure,
    subplot(1,3,1)
    hold on
    plot((params.r+5)*sin(t),(params.r+5)*cos(t),'k', 'linewidth', 2)
    plot(dt.X,dt.Y)
    axis square
    axis([-5-params.r 5+params.r -5-params.r 5+params.r])
    subplot(1,3,[2 3])
    hold on
    plot((1:(length(dt.Vr)-1))/60, zeros((length(dt.Vr)-1),1), 'color', [0.6 0.6 0.6])
    plot((1:(length(dt.Vr)-1))/60,dt.Vr(1:end-1), 'color', [0 0 1])
    plot((1:(length(dt.Vf)-1))/60,-1000+10*dt.Vt(1:end-1), 'k')
    plot((1:(length(fM)))/60, 1000+60*(fM(1:end)), 'color', [1 0 0])
    plot((1:(length(fV)))/60, 1000+60*(fV(1:end)), 'color', [0.5 0 0])
    plot((1:(length(PER)))/60, 1000-300*abs(PER(1:end)), 'color', [0 0 0.5])
    for i = 1 : length(dt.Bouts)
        plot(dt.Bouts{i}/60, dt.Vr(dt.Bouts{i}), 'g')
        plot(dt.Bouts{i}/60,-1000+10*dt.Vt(dt.Bouts{i}), 'g')
    end
    for i = 1 : length(str.FBouts)
        plot(str.FBouts{i}/60, dt.Vr(str.FBouts{i}), 'r')
        plot(str.FBouts{i}/60,-1000+10*dt.Vt(str.FBouts{i}), 'r')
    end
    axis([0 30 -1300 1800])
end

if params.dispStr == 1
    disp('-----------------------------------------------')
    disp(['str = ' num2str(str.STR'*str.N/sum(str.N))])
    disp(['vrBTF = ' num2str(pss.vrl'*pss.nvrl/sum(pss.nvrl))])
    disp(['vrFTB = ' num2str(pss.vrr'*pss.nvrr/sum(pss.nvrr))])
    disp(['pssBTF = ' num2str(pss.pvrl'*pss.nvrl/sum(pss.nvrl))])
    disp(['pssFTB = ' num2str(pss.pvrr'*pss.nvrr/sum(pss.nvrr))])
    disp(['nBTF = ' num2str(sum(pss.nvrl))])
    disp(['nFTB = ' num2str(sum(pss.nvrr))])
    disp(['dst = ' num2str(str.DST'*str.N/sum(str.N))])
    disp(['pss = ' num2str(pss.pss)])
    disp('-----------------------------------------------')
end
end