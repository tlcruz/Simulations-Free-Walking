function [sim_data] = Andre_2QuadModel(eye_sample, ...
    lp_Tau_HR, hp_Tau_HR,fs_LED,rec,hwr,rightWeight,onWeight)
% simulate the flight arena, requires eye_filt map, the Pattern to display,
% and the time series that specifies the frame positions. Also need to know
% the sample rate (as fps). Can specify a period of blank display, by
% setting values in frame_positions to -1, during these period, display
% will show intermediate value (no apparent motion). Also send in tc in
% seconds.
% this version now runs 2 half-eye EMD, to make sure all is symmetric

% [~, num_receptors] = size(eye_sample);
% % how many receptors per eye?
% rec_pe = num_receptors/2; % currently assume same number per eye,
% % deal with separately if this is not the case

% initializations for HR model

%It's constants time!!!!
%general stuff
DC = 10;%10;
% hwr = 0.05;
dt = 1/fs_LED;

% for weigth4 = -1:0.1:1
%weights
weigth1 = rightWeight; % right detector
weigth2 = 1-rightWeight; % left detector
weigth3 = onWeight; % ON pathway
weigth4 = 1-onWeight; % OFF pathway

TfiltLP = dt/lp_Tau_HR;
TfiltHP = dt/hp_Tau_HR;
h = 1/fs_LED;  % the sampling interval
% A_lp = 1 - (2*lp_Tau_HR)/h; B_lp = 1 + (2*lp_Tau_HR)/h;  % the 2 filter coefficients
% %%Here we use a bilinear transform

% A_hp = hp_Tau_HR/(hp_Tau_HR+h);
%
% HR_Motion = zeros(num_frames, num_receptors - 2);   %% due to motion detectors at the end


% InMat     = 5*(rand(1,num_receptors) - 0.5);    % input into eye
% InMat_1   = 5*(rand(1,num_receptors) - 0.5);    % last input value for causal filter
% FiltMat   = zeros(size(InMat));
% FiltMat_1 = zeros(size(InMat));                 % filter
% these start off with random numbers

%lamina cells
in = eye_sample;
hp = filter([1-TfiltHP TfiltHP-1],[1 TfiltHP-1],in);
dc = in*DC/100;
f = hp+dc;
% f = in;
%half-wave rectifier
% ON rectifier, with zero
l1 = f.*gt(f,0);


%OFF rectifier, where hwr is taken into account. See Eichner, Borst 2011
l2 = f.*lt(f,0);


num_receptors = rec;

if rec == 61 || rec == 72;
    for j = 1:rec-1
        
        % for j = 4:4; %use this for just 2 photoreceptors
        % for j = 4:5; %use this for just 3 photoreceptors
        
        
        % EMD model
        %ON pathway
        % for tauLP = 0.01:0.01:1
        % FOR FILTERING
        % http://stackoverflow.com/questions/1783633/matlab-apply-aLP-low-pass-or-high-pass-filter-to-an-array
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
%         HR_Motion(:,j) = weigth3*dOFF(j,:);
        
        %     response2(i,:) = sim('Simulink_attempt');
    end

elseif rec == 10
    for j = 1:10
        
        % for j = 4:4; %use this for just 2 photoreceptors
        % for j = 4:5; %use this for just 3 photoreceptors
        
        
        % EMD model
        %ON pathway
        % for tauLP = 0.01:0.01:1
        % FOR FILTERING
        % http://stackoverflow.com/questions/1783633/matlab-apply-aLP-low-pass-or-high-pass-filter-to-an-array
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
        
        %     response2(i,:) = sim('Simulink_attempt');
    end
    
    
elseif rec == 3
    for j = 4:5; %use this for just 3 photoreceptors
        
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
        
        %     response2(i,:) = sim('Simulink_attempt');
    end
    
    
elseif rec == 2
    for j = 4:4; %use this for just 2 photoreceptors
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
        %     response2(i,:) = sim('Simulink_attempt');
    end
end





%
% figure
% % subplot(2,1,1)
% plot(t_LED(1:10*fs_LED), sim_data.HR_mean_ts(1:10*fs_LED));
% title('Response with the 2 quadrant detector (One Stripe - Right)')
% ylabel('Membrane Potential (arbitraty unit)')
% set(gca,'FontSize', 16)
%
% end
%     Tuthill
%     InMat = A_hp*(FiltMat_1)+ A_hp*(InMate-InMat_1);
%     %%y(n-1) is the previous filtered output
%     %%x(n) is the current, unfiltered input
%      %%x(n-1) is the previous filtered input
%
%     FiltMat = ( InMate + InMat_1 - A_lp*FiltMat_1 ) / B_lp ;
%     %%y(n-1) is the previous, lp filtered input (FiltMat_1)
%     %%x(n) is the current unfiltered input  (Inmate)
%      %%x(n-1) is the previous filtered input (Inmat_1)
%
%     %%Add signal and previous signal and subtract filtered previous input;
%
%     InMat_1   = InMat; FiltMat_1 = f1iltMat;      %%resets these for next round
%     HR_Motion(j,1:(rec_pe-1)) = (FiltMat(1:(rec_pe-1)).*InMat(2:rec_pe) - FiltMat(2:rec_pe).*InMat(1:(rec_pe-1)));    %%correlate and subtracts reichardt thing
%     HR_Motion(j,(rec_pe):(2*rec_pe - 2)) = -((FiltMat((rec_pe+2):end).*InMat((rec_pe+1):end-1) - FiltMat((rec_pe+1):end-1).*InMat((rec_pe+2):end)));
%     %HR_Motion(j,:) = (FiltMat(1:end-1).*InMat(2:end) - FiltMat(2:end).*InMat(1:end-1));

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
