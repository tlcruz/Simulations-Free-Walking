function PlotSingleSim(outp, ti, tf)
if nargin < 3
   ti = 0;
   tf = (length(outp.dt.Vt)-1)/60;
end
params = outp.params;
PER = outp.PER;
pertVrMat = outp.pertVrMat;
pertMat = outp.pertMat;
strA = outp.strA;
str = outp.str;
angd = outp.angd;
dist = outp.dist;
N = outp.n;
pss = outp.pss;
npss = outp.npss;
vrl = outp.vrl;
pvrl = outp.pvrl;
nvrl = outp.nvrl;
vrr = outp.vrr;
pvrr = outp.pvrr;
nvrr = outp.nvrr;
dt = outp.dt;
fM = outp.fM;
fV = outp.fV;
t = -0.1:0.1:2*pi;

fi = ti*60+1;
ff = tf*60;
figure,
subplot(1,3,1)
hold on
plot((params.r+5)*sin(t),(params.r+5)*cos(t),'k', 'linewidth', 2)
plot(dt.X(fi:ff),dt.Y(fi:ff))
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
for i = 1 : length(strA.FBouts)
    plot(strA.FBouts{i}/60, dt.Vr(strA.FBouts{i}), 'r')
    plot(strA.FBouts{i}/60,-1000+10*dt.Vt(strA.FBouts{i}), 'r')
end
axis([ti tf -1300 1800])

disp('-----------------------------------------------')
disp(['str = ' num2str(str)])
disp(['vrBTF = ' num2str(vrl'*nvrl/sum(nvrl))])
disp(['vrFTB = ' num2str(vrr'*nvrr/sum(nvrr))])
disp(['pssBTF = ' num2str(pvrl'*nvrl/sum(nvrl))])
disp(['pssFTB = ' num2str(pvrr'*nvrr/sum(nvrr))])
disp(['nBTF = ' num2str(sum(nvrl))])
disp(['nFTB = ' num2str(sum(nvrr))])
disp(['dst = ' num2str(dist)])
disp(['pss = ' num2str(pss)])
disp('-----------------------------------------------')