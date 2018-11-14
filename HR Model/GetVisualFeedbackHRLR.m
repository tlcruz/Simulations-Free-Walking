function [vF, simLeft, simRight, params] = GetVisualFeedbackHRLR(Vr, params)
% Visual Stimulus
inputArray = zeros(params.inputSize,1);
inputArray(1:ceil(params.D*params.inputSize)) = 1;
xi = 1 : params.inputSize;
auxArray = downsample(inputArray,params.AngSize);
auxxi = downsample(xi, params.AngSize);
vFL = [];
vFR = [];
for k = 1 : params.NCopies
    auxArray = auxArray(randperm(length(auxArray)));
    if params.AngSize > 4
        while ~isempty(find(diff(find(auxArray==1))==1))
            auxArray = auxArray(randperm(length(auxArray)));
        end
    end
    %     RandInput = vertcat(RandInput, interp1(auxxi, auxArray, xi, 'nearest'));
    RandInput = interp1(auxxi, auxArray, xi, 'nearest');
    RandInput(isnan(RandInput)) = 0 ;
    
    % Input Velocity
    tVel = length(Vr)/params.fsVel;
    inputVel = round((180/pi)*Vr/params.fsLED);
    in = zeros(floor(tVel*params.fsLED), params.inputSize);
    
    % Generate Full Range Stimulus
    in(1,:) = RandInput(1,:);
    vv = [];
    for k = 2: size(in,1)
        in(k,:) = circshift(in(k-1,:)',...
            [round(inputVel(round(k*params.fsVel/params.fsLED))) 0]);
        vv = vertcat(vv, round(inputVel(round(k*params.fsVel/params.fsLED))));
    end
    
    % Smooth with photoreceptor filter
    resized_in = imresize(in, [tVel*params.fsLED 72]);
    
    % Divide in Left and Right eye
    inputLeft = (1-params.WeightL)*resized_in(:, 1:floor(end/2));
    inputRight = resized_in(:, 1+floor(end/2):end);
    
    % Apply the model
    [simLeft] = DQuadModel(inputLeft,15e-3,50e-3,60,...
        36,0,1/2,1/2);
    [simRight] = DQuadModel(inputRight,15e-3,50e-3,60,...
        36,0,1/2,1/2);
    
    vFL = vertcat(vFL, mean(simLeft.HR_mean_ts(end-3:end)));
    vFR = vertcat(vFR, mean(simRight.HR_mean_ts(end-3:end)));
end
vFL = 1000*mean(vFL);
vFR = 1000*mean(vFR);

vFFTBL = 0;
vFFTBR = 0;
vFBTFL = 0;
vFBTFR = 0;

% params.WeightBTF = params.ab/params.vf;
% params.WeightFTB = params.vf/(20+params.ab);

params.WeightBTF = 0.6;
params.WeightFTB = 0.6;

if vFL > 0
   vFBTFL = vFL*params.WeightBTF;
   vFL = vFL*params.WeightBTF;
else
   vFFTBL = vFL*params.WeightFTB; 
   vFL = vFL*params.WeightFTB; 
end
if vFR > 0
   vFFTBR = vFR*params.WeightFTB;
   vFR = vFR*params.WeightFTB;
else
   vFBTFR = vFR*params.WeightBTF; 
   vFR = vFR*params.WeightBTF; 
end

vFBTF = vFBTFL+vFBTFR;
vFFTB = vFFTBL+vFFTBR;

vF.vFFTBL = vFFTBL;
vF.vFFTBR = vFFTBR;
vF.vFFTB = vFFTB;
vF.vFBTFL = vFBTFL;
vF.vFBTFR = vFBTFR;
vF.vFBTF = vFBTF;
vF.vFL = vFL;
vF.vFR = vFR;
end


















