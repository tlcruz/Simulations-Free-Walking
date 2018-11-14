function [vF, sim] = GetVisualFeedbackHR(Vr, params)
% Visual Stimulus
inputArray = zeros(params.inputSize,1);
inputArray(1:ceil(params.D*params.inputSize)) = 1;
xi = 1 : params.inputSize;
auxArray = downsample(inputArray,params.AngSize);
auxxi = downsample(xi, params.AngSize);
vF = [];
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
        vv = vertcat(vv, ...
            round(inputVel(round(k*params.fsVel/params.fsLED))));
    end
    
    % Smooth with photoreceptor filter
    resized_in = imresize(in, [tVel*params.fsLED 72]);
    
    % Apply the model
    [sim] = DQuadModel(resized_in,15e-3,50e-3,60,...
        72,0,1/2,1/2);
    
    vF = vertcat(vF, mean(sim.HR_mean_ts(end-3:end)));
end
vF = 1000*mean(vF);
end


















