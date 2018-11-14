

%% Used Angular Sizes: 1, 3, 5, 10, 15

load('Dt.mat')


k =1;

fsLED = 60; % LED frequency. 
num_receptors = 72; % number of receptors

D = 0.2;
inputSize = 270;
NCopies = 10;
Vel = Dt.VrNG{1,1}{k,1};
fsVel = 60;
tVel = length(Vel)/fsVel;
AngSize = 1;

inputVel = zeros(1,tVel*fsLED);
inputVel = round(Vel/fsLED);

inputArray = zeros(inputSize,1);
inputArray(1:ceil(D*inputSize)) = 1;
xi = 1 : inputSize;
auxArray = downsample(inputArray,AngSize);
auxxi = downsample(xi, AngSize);
RandInput = [];
for k = 1 : NCopies
    auxArray = auxArray(randperm(length(auxArray)));
    RandInput = vertcat(RandInput, interp1(auxxi, auxArray, xi, 'nearest'));
end
figure,
imagesc(RandInput)

RandInput(isnan(RandInput)) = 0 ;
in = zeros(tVel*fsLED,inputSize);
in(1,:) = RandInput(1,:);
%%
vv = [];
for k = 2: length(in)
    in(k,:) = circshift(in(k-1,:)',[round(inputVel(round(k*fsVel/fsLED))) 0]);
    vv = vertcat(vv, round(inputVel(round(k*fsVel/fsLED))));
end

figure; imagesc(in)

%% Smooth with photoreceptor filter.


resized_in = imresize(in, [tVel*fsLED 72 ]);


figure;
subplot(2,1,1)
imagesc(in)
subplot(2,1,2)
imagesc(resized_in)


%% Apply the model
%                 Original from Tuthill
lp_Tau_HR = 15e-3;  % time constant of the lp-filter
hp_Tau_HR = 50e-3;  % time constant of the hp filter, from Borst et al, 2003

hwr = 0.0;
rightWeight = 1/2;
onWeight = 1/2;
rec = 72;

% Andre_make_eye_filters; % currently this is w. 4.6 deg ommatidia, 72 total

[resized_sim_data2Quad] = Andre_2QuadModel ...
    (resized_in, lp_Tau_HR, hp_Tau_HR,fsLED,rec,hwr,rightWeight,onWeight)';

[sim_data2Quad] = Andre_2QuadModel ...
    (in, lp_Tau_HR, hp_Tau_HR,fsLED,rec,hwr,rightWeight,onWeight)';
%%
figure; 
% plot(sim_data2Quad.HR_mean_ts)
% hold on
plot(resized_sim_data2Quad.HR_mean_ts)
hold on
plot(0.001*vv)
legend('resized','Vr')
% legend('resized','Vr')

