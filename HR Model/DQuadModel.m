function [sim_data] = DQuadModel(eye_sample, ...
    lp_Tau_HR, hp_Tau_HR,fs_LED,rec,hwr,rightWeight,onWeight)

if nargin == 1
    lp_Tau_HR = 15e-3;
    hp_Tau_HR = 50e-3;
    fs_LED = 60;
    rec = 72;
    hwr = 0.0;
    rightWeight = 1/2;
    onWeight = 1/2;
end

DC = 10;
dt = 1/fs_LED;
TfiltLP = dt/lp_Tau_HR;
TfiltHP = dt/hp_Tau_HR;

%weights
weigth1 = rightWeight; % right detector
weigth2 = 1-rightWeight; % left detector
weigth3 = onWeight; % ON pathway
weigth4 = 1-onWeight; % OFF pathway

%lamina cells
in = eye_sample;
hp = filter([1-TfiltHP TfiltHP-1],[1 TfiltHP-1],in);
dc = in*DC/100;
f = hp+dc;

%half-wave rectifier
% ON rectifier, with zero
l1 = f.*gt(f,0);
%OFF rectifier, where hwr is taken into account. See Eichner, Borst 2011
l2 = f.*lt(f,0);

if rec > 3;
    for j = 1:rec-1
        lpON_1 = filter(TfiltLP,[1 TfiltLP-1],l1(:,j));
        lpON_2 = filter(TfiltLP,[1 TfiltLP-1],l1(:,j+1));
        
        %Multiplication
        dON_1 = lpON_1 .* l1(:,j+1);
        dON_2 = l1(:,j) .* lpON_2;
        
        %Subtraction
        dON(j,:) = weigth1*dON_1 - weigth2*dON_2;
        
        %OFF pathway
        lpOFF_1 = filter(TfiltLP,[1 TfiltLP-1],l2(:,j));
        lpOFF_2 = filter(TfiltLP,[1 TfiltLP-1],l2(:,j+1));
        
        
        %Multiplication
        dOFF_1 = lpOFF_1 .* l2(:,j+1);
        dOFF_2 = l2(:,j) .* lpOFF_2;
        
        %Subtraction
        dOFF(j,:) = weigth1*dOFF_1 - weigth2*dOFF_2;
        
        % Sum
        HR_Motion(:,j) = weigth3*dON(j,:) + weigth4*dOFF(j,:);

    end
end

sim_data.eye_sample = eye_sample;
sim_data.HR_Motion = HR_Motion;
% save the results for HS and different things
sim_data.HR_mean_ss = mean(sim_data.HR_Motion);
sim_data.HR_mean_ss_avg =  mean(sim_data.HR_mean_ss);
sim_data.HR_mean_ts = mean(sim_data.HR_Motion, 2);
sim_data.HR_sum_ts = sum(sim_data.HR_Motion, 2);
sim_data.HR_max_ts = max(sim_data.HR_mean_ts);
sim_data.HR_mean_eye = mean(sim_data.eye_sample, 2);

% save the results for T4T5
sim_data.T4T5_dON = dON;
sim_data.T4T5_dOFF = dOFF;
end
