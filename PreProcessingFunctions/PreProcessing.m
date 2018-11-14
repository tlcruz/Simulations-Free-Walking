function [dt] = PreProcessing(XX,YY,TH)
addpath('PreProcessingFunctions');

VfRaw = 60*((XX(2:end)-XX(1:end-1)).*cos(pi*TH(1:end-1)/180) + ...
    (YY(2:end)-YY(1:end-1)).*sin(pi*TH(1:end-1)/180));
VsRaw = 60*(-(XX(2:end)-XX(1:end-1)).*sin(pi*TH(1:end-1)/180) + ...
    (YY(2:end)-YY(1:end-1)).*cos(pi*TH(1:end-1)/180));
VrRaw = 60*diff(TH);
Angle = TH(1:end);
WallDist = 45 - sqrt(XX.*XX + YY.*YY);
[VfRaw, VsRaw, VrRaw, VtRaw, Angle, nJ] = ClearJumps(VfRaw, VsRaw, VrRaw, Angle);
if nJ ~= 0
   disp('Error nJ') 
end

nPoints = 6;
Vr = smooth(vertcat(VrRaw,0), nPoints/length(VrRaw), 'lowess');
Vf = smooth(vertcat(VfRaw,0), nPoints/length(VfRaw), 'lowess');
Vt = smooth(vertcat(VtRaw,0), nPoints/length(VtRaw), 'lowess');
Vs = smooth(vertcat(VsRaw,0), nPoints/length(VsRaw), 'lowess');

[actState, ActBouts] = GetBouts(Vr, Vf, Vs);


dt.Vrus = VrRaw;
dt.Vf = Vf;
dt.Vs = Vs;
dt.Vr = Vr;
dt.Vt = Vt;
dt.Ang = Angle;
dt.X = XX;
dt.Y = YY;
dt.WD = WallDist;
dt.actSt = actState;
dt.Bouts = ActBouts;

end