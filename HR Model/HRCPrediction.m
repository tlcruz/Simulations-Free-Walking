function [Vresp, Vrp, Vmean, Vstd] = HRCPrediction(Vr, ang, NCopies)
fsLED = 60; % LED frequency. 
num_receptors = 72; % number of receptors
D = 0.20;%0.17;
inputSize = 360;
fsVel = 60;
lp_Tau_HR = 15e-3;%15e-3;  % time constant of the lp-filter
hp_Tau_HR = 50e-3;%50e-3;  % time constant of the hp filter, from Borst et al, 2003
hwr = 0.0;
rightWeight = 1/2;
onWeight = 1/2;
rec = 72;

tVel = length(Vr)/fsVel;
inputVr = round(Vr/fsLED);

inputArray = zeros(inputSize,1);
inputArray(1:ceil(D*inputSize)) = 1;
xi = 1 : inputSize;
auxArray = downsample(inputArray,ang);
auxxi = downsample(xi, ang);
RandInput = [];
for k = 1 : NCopies
    auxArray = auxArray(randperm(length(auxArray)));
    RandInput = vertcat(RandInput, interp1(auxxi, auxArray, xi, 'nearest'));
end
RandInput(isnan(RandInput)) = 0 ;
in = zeros(floor(tVel*fsLED),inputSize);

Vresp = [];
for n = 1 : NCopies
in(1,:) = RandInput(n,:);
vv = [];
for k = 2: size(in,1)
    in(k,:) = circshift(in(k-1,:)',[round(inputVr(ceil(k*fsVel/fsLED))) 0]);
    vv = vertcat(vv, round(inputVr(ceil(k*fsVel/fsLED))));
end
Vrp = vv;
resized_in = imresize(in, [floor(tVel*fsLED) 72 ]);
[resized_sim_data2Quad] = Andre_2QuadModel ...
    (resized_in, lp_Tau_HR, hp_Tau_HR,fsLED,rec,hwr,rightWeight,onWeight)';
Vresp = horzcat(Vresp,resized_sim_data2Quad.HR_mean_ts);
end

Vmean = mean(Vresp,2);
Vstd = std(Vresp,1,2);
end