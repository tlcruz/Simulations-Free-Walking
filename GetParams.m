function [params] = GetParams()
    % display params 
    params.plotTraj = 0;
    params.dispStr = 0;

    % Trajectory Simulation Parameters
    params.x0 = 0;
    params.y0 = 0;
    params.th0 = 0;
    params.L = 3*60*60;
    params.vf = 20;
    params.noisesize = 20;
    params.r = 35;
    fc = 15;
    fs = 60;
    [params.noisebf,params.noiseaf] = butter(6,fc/(fs/2));
%     params.spkPa = 0.002;
    params.spkPa = 0.0083;
    params.spkPb = 0.001;
    params.spkPl = 0.16;
    params.spkAa = 200;
    params.spkAb = 2000;
    params.spkAl = 0.1872/3;
    params.spkAdx = 50;
    
    % Control
    params.pathKernelMotor = 'KernelMotorFeedback.mat';
    aux = load(params.pathKernelMotor);
    params.kernelMotor = aux.f;
    params.kernelMotorSize = 5;

    % params initial peak detection
    params.Fs = 60;
    params.minF = 10;
    params.maxF = 15;
    params.nstds = 2; % 
    params.MinPeakDist = 0;
    params.thrDist = 5;
    params.MinPeakPromF = 0;
    params.MinPeakPromV = 20;
    params.maxP = 9;
    params.minP = 7;
    
    % params FBuots
    params.pkvt = 200;
    params.promt = 150;
    params.vfstdt = 3;
    params.vfmt = 6;
    params.vft = 6;
    
    % params Straightness
    params.spkTempPath = 'SpikeTemplateL.mat';
    params.cutoff = 0.20;
    params.RG = 0;
    params.thr = 2*60;
    params.thrspk = 100;
    params.delta = 15;
    params.sep = 20;
    params.windStr = 20;
    params.maxvf = 0.4;
    params.minvf = 5;
    params.minStrB = 20; %
    params.mDistWall = 0;
    params.WindowCCVFW = 20;
    params.WindowCCVALL = 40;
    params.btthr = 21;
    
    % params HR Model
    params.fsLED = 60;
    params.D = 0.2;
    params.inputSize = 360;
    params.NCopies = 1;
    params.AngSize = 10;
    params.fsVel = 60;
    params.lp_Tau_HR = 15e-3;  % time constant of the lp-filter
    params.hp_Tau_HR = 50e-3;  % time constant of the hp filter, from Borst et al, 2003
    params.hwr = 0.0;
    params.rightWeight = 1/2;
    params.onWeight = 1/2;
    params.rec = 72;
    params.kVisualSize = 10;
    params.WeightR = 1;
    params.WeightL = 1;
    
    
    params.ab = 5;
    
    % perturbations
    params.WeightBTF = 1;
    params.WeightFTB = 1;
    params.PerP = 'FTBR';
    params.SizePerturbation = 4;
    params.PertOn = 1;
    params.pertD = 10;
end