[params] = GetParams();
params.noisesize = 15;
params.kernelMotorSize = 4.9;
params.NCopies = 1;
params.kVisualSize = -6;
params.AngSize = 10;
params.delta = 15;
params.L = 1*60*60;
params.ab = 8.5;
% params.WeightBTF = 0.3;
% params.WeightFTB = 0.5;
params.PertOn = 0;
outp = StartSimulation(params);
PlotSingleSim(outp, 9, 30);

%%
for dotSize = [1 2 3 5 10]
    poolobj = gcp('nocreate');
    % ds = [10];%[1 2 3 5 10];
    % ds = [0 0.25 0.5 0.75 1];
    ds = 0.25:0.25:1;
    GMS = cell(length(ds), 1);
    STDS = cell(length(ds), 1);
    GMPSS = cell(length(ds), 2);
    STDPSS = cell(length(ds), 2);
    GMADL = cell(length(ds), 2);
    STDADL = cell(length(ds), 2);
    GMADR = cell(length(ds), 2);
    STDADR = cell(length(ds), 2);
    GMPL = cell(length(ds), 2);
    STDPL = cell(length(ds), 2);
    GMPR = cell(length(ds), 2);
    STDPR = cell(length(ds), 2);
    if isempty(poolobj)
        parpool(12)
    end
    nC = 12;
    a = 1;
    for i = ds
        GMM = cell(12,1);
        GN = cell(12,1);
        PSS1M = cell(12,1);
        PSS1N = cell(12,1);
        PSS2M = cell(12,1);
        PSS2N = cell(12,1);
        AL1 = cell(12,1);
        AL2 = cell(12,1);
        PL1 = cell(12,1);
        PL2 = cell(12,1);
        NAL1 = cell(12,1);
        NAL2 = cell(12,1);
        AR1 = cell(12,1);
        AR2 = cell(12,1);
        PR1 = cell(12,1);
        PR2 = cell(12,1);
        NAR1 = cell(12,1);
        NAR2 = cell(12,1);
        parfor k = 1 : 12
            [params] = GetParams();
            params.kernelMotorSize = 5;
            params.noisesize = 20;
            params.kVisualSize = 12;
            params.AngSize = dotSize;
            params.WeightR = i;
            params.NCopies = 10;
            [outp1] = StartSimulation(params);
            params.kVisualSize = -12;
            [outp2] = StartSimulation(params);
            GMM{k} = outp1.str;
            GN{k} = outp1.n;
            AL1{k} = outp1.vrl;
            PL1{k} = outp1.pvrl;
            NAL1{k} = outp1.nvrl;
            AL2{k} = outp2.vrl;
            PL2{k} = outp2.pvrl;
            NAL2{k} = outp2.nvrl;
            AR1{k} = outp1.vrr;
            PR1{k} = outp1.pvrr;
            NAR1{k} = outp1.nvrr;
            AR2{k} = outp2.vrr;
            PR2{k} = outp2.pvrr;
            NAR2{k} = outp2.nvrr;
            PSS1M{k} = outp1.pss;
            PSS1N{k} = outp1.npss;
            PSS2M{k} = outp2.pss;
            PSS2N{k} = outp2.npss;
        end
        GMM = cell2mat(GMM);
        GN = cell2mat(GN);
        gm = GMM'*GN/sum(GN);
        sem = sqrt(((GMM-gm).*(GMM-gm))'*GN/sum(GN))/sqrt(10);
        GMS{a} = gm;
        STDS{a} = sem;
        GMPSS1M = cell2mat(PSS1M);
        GMPSS1N = cell2mat(PSS1N);
        gmpss1 = GMPSS1M'*GMPSS1N/sum(GMPSS1N);
        sempss1 = sqrt(((GMPSS1M-gmpss1).*(GMPSS1M-gmpss1))'*GMPSS1N/sum(GMPSS1N))/sqrt(10);
        GMPSS{a,1} = gmpss1;
        STDPSS{a,1} = sempss1;
        AL1 = cell2mat(AL1);
        NAL1 = cell2mat(NAL1);
        gmal1 = AL1'*NAL1/sum(NAL1);
        semal1 = sqrt(((AL1-gmal1).*(AL1-gmal1))'*NAL1/sum(NAL1))/sqrt(nC);
        GMADL{a,1} = gmal1;
        STDADL{a,1} = semal1;
        
        PL1 = cell2mat(PL1);
        gmpl1 = PL1'*NAL1/sum(NAL1);
        sempl1 = sqrt(((PL1-gmpl1).*(PL1-gmpl1))'*NAL1/sum(NAL1))/sqrt(nC);
        GMPL{a,1} = gmpl1;
        STDPL{a,1} = sempl1;
        
        AR1 = cell2mat(AR1);
        NAR1 = cell2mat(NAR1);
        gmar1 = AR1'*NAR1/sum(NAR1);
        semar1 = sqrt(((AR1-gmar1).*(AR1-gmar1))'*NAR1/sum(NAR1))/sqrt(nC);
        GMADR{a,1} = gmar1;
        STDADR{a,1} = semar1;
        
        PR1 = cell2mat(PR1);
        gmpr1 = PR1'*NAR1/sum(NAR1);
        sempr1 = sqrt(((PR1-gmpr1).*(PR1-gmpr1))'*NAR1/sum(NAR1))/sqrt(nC);
        GMPR{a,1} = gmpr1;
        STDPR{a,1} = sempr1;
        
        AL2 = cell2mat(AL2);
        NAL2 = cell2mat(NAL2);
        gmal2 = AL2'*NAL2/sum(NAL2);
        semal2 = sqrt(((AL2-gmal2).*(AL2-gmal2))'*NAL2/sum(NAL2))/sqrt(nC);
        GMADL{a,2} = gmal2;
        STDADL{a,2} = semal2;
        
        PL2 = cell2mat(PL2);
        gmpl2 = PL2'*NAL2/sum(NAL2);
        sempl2 = sqrt(((PL2-gmpl2).*(PL2-gmpl2))'*NAL2/sum(NAL2))/sqrt(nC);
        GMPL{a,2} = gmpl2;
        STDPL{a,2} = sempl2;
        
        AR2 = cell2mat(AR2);
        NAR2 = cell2mat(NAR2);
        gmar2 = AR2'*NAR2/sum(NAR2);
        semar2 = sqrt(((AR2-gmar2).*(AR2-gmar2))'*NAR2/sum(NAR2))/sqrt(nC);
        GMADR{a,2} = gmar2;
        STDADR{a,2} = semar2;
        
        PR2 = cell2mat(PR2);
        gmpr2 = PR2'*NAR2/sum(NAR2);
        sempr2 = sqrt(((PR2-gmpr2).*(PR2-gmpr2))'*NAR2/sum(NAR2))/sqrt(nC);
        GMPR{a,2} = gmpr2;
        STDPR{a,2} = sempr2;
        
        GMPSS2M = cell2mat(PSS2M);
        GMPSS2N = cell2mat(PSS2N);
        gmpss2 = GMPSS2M'*GMPSS2N/sum(GMPSS2N);
        sempss2 = sqrt(((GMPSS2M-gmpss2).*(GMPSS2M-gmpss2))'*GMPSS2N/sum(GMPSS2N))/sqrt(10);
        GMPSS{a,2} = gmpss2;
        STDPSS{a,2} = sempss2;
        disp(['Dot Size ' num2str(i) ' Str: ' num2str(gm) ...
            ' +/- ' num2str(sem) '  VI:' num2str((gmpss2-gmpss1)/(gmpss2+gmpss1))])
        a = a + 1;
    end
    dtt.GMS = GMS;
    dtt.STDS = STDS;
    dtt.GMPSS = GMPSS;
    dtt.STDPSS = STDPSS;
    dtt.GMADL = GMADL;
    dtt.STDADL = STDADL;
    dtt.GMADR = GMADR;
    dtt.STDADR = STDADR;
    dtt.GMPL = GMPL;
    dtt.STDPL = STDPL;
    dtt.GMPR = GMPR;
    dtt.STDPR = STDPR;
    save(['C:\Users\tomas\Desktop\dt2nm' num2str(dotSize) '.mat'], 'dtt');
end
%%
load('C:\Users\tomas\Desktop\dt210.mat')
ds = 0:0.1:1;
figure,
cmap = autumn(length(ds));
cmap2 = winter(length(ds));
for i = 1:length(ds)
    subplot(1,5,1)
    errorbar(ds(i), dtt.GMS{i}, dtt.STDS{i}, 'color', cmap(i,:))
    hold on
    ylabel('Straightness')
    axis([(min(ds)-0.1) (max(ds)+0.1) 20 60])
    
    subplot(1,5,2)
    errorbar(ds(i), dtt.GMPL{i,1}, dtt.STDPL{i,1}, 'color', cmap2(i,:))
    hold on
    errorbar(ds(i), dtt.GMPR{i,2}, dtt.STDPR{i,2}, 'color', cmap(i,:))
    ylabel('FTB?')
    axis([(min(ds)-0.1) (max(ds)+0.1) 0.35 0.95])
    
    subplot(1,5,3)
    errorbar(ds(i), dtt.GMPR{i,1}, dtt.STDPR{i,1}, 'color', cmap2(i,:))
    hold on
    errorbar(ds(i), dtt.GMPL{i,2}, dtt.STDPL{i,2}, 'color', cmap(i,:))
    ylabel('BTF?')
    axis([(min(ds)-0.1) (max(ds)+0.1) 0.35 0.95])
    
    subplot(1,5,4)
    errorbar(ds(i), dtt.GMADL{i,1}, dtt.STDADL{i,1}, 'color', cmap2(i,:))
    hold on
    errorbar(ds(i), dtt.GMADR{i,2}, dtt.STDADR{i,2}, 'color', cmap(i,:))
    ylabel('FTB?')
    axis([(min(ds)-0.1) (max(ds)+0.1) -500 700])
    
    subplot(1,5,5)
    errorbar(ds(i), dtt.GMADR{i,1}, dtt.STDADR{i,1}, 'color', cmap2(i,:))
    hold on
    errorbar(ds(i), dtt.GMADL{i,2}, dtt.STDADL{i,2}, 'color', cmap(i,:))
    ylabel('BTF?')
    axis([(min(ds)-0.1) (max(ds)+0.1) -500 700])
end
%%
d1 = load('C:\Users\tomas\Desktop\dt2nm1.mat');
d1 = d1.dtt;
d2 = load('C:\Users\tomas\Desktop\dt2nm2.mat');
d2 = d2.dtt;
d3 = load('C:\Users\tomas\Desktop\dt2nm3.mat');
d3 = d3.dtt;
d5 = load('C:\Users\tomas\Desktop\dt2nm5.mat');
d5 = d5.dtt;
d10 = load('C:\Users\tomas\Desktop\dt2nm10.mat');
d10 = d10.dtt;
dotSize = [1 2 3 5 10];
i = 2;
figure,
subplot(1,2,1)
hold on
errorbar(dotSize(1), d1.GMPL{i,1}, d1.STDPL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(1), d1.GMPR{i,2}, d1.STDPR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
errorbar(dotSize(2), d2.GMPL{i,1}, d2.STDPL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(2), d2.GMPR{i,2}, d2.STDPR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
errorbar(dotSize(3), d3.GMPL{i,1}, d3.STDPL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(3), d3.GMPR{i,2}, d3.STDPR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
errorbar(dotSize(4), d5.GMPL{i,1}, d5.STDPL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(4), d5.GMPR{i,2}, d5.STDPR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
errorbar(dotSize(5), d10.GMPL{i,1}, d10.STDPL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(5), d10.GMPR{i,2}, d10.STDPR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
hold on
errorbar(dotSize(1), d1.GMPR{i,1}, d1.STDPR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(1), d1.GMPL{i,2}, d1.STDPL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
errorbar(dotSize(2), d2.GMPR{i,1}, d2.STDPR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(2), d2.GMPL{i,2}, d2.STDPL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
errorbar(dotSize(3), d3.GMPR{i,1}, d3.STDPR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(3), d3.GMPL{i,2}, d3.STDPL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
errorbar(dotSize(4), d5.GMPR{i,1}, d5.STDPR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(4), d5.GMPL{i,2}, d5.STDPL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
errorbar(dotSize(5), d10.GMPR{i,1}, d10.STDPR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(5), d10.GMPL{i,2}, d10.STDPL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
ylabel('PSS FTB - BTF')
xlabel('Dot Size')
set(gca,'xtick',[1 2 3 5 10]);
axis([0 11 0.4 0.9])
subplot(1,2,2)
hold on
errorbar(dotSize(1), d1.GMADL{i,1}, d1.STDADL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(1), d1.GMADR{i,2}, d1.STDADR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
errorbar(dotSize(2), d2.GMADL{i,1}, d2.STDADL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(2), d2.GMADR{i,2}, d2.STDADR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
errorbar(dotSize(3), d3.GMADL{i,1}, d3.STDADL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(3), d3.GMADR{i,2}, d3.STDADR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
errorbar(dotSize(4), d5.GMADL{i,1}, d5.STDADL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(4), d5.GMADR{i,2}, d5.STDADR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
errorbar(dotSize(5), d10.GMADL{i,1}, d10.STDADL{i,1},'o', 'color', 'r','MarkerSize',10, ...
    'MarkerFaceColor','r')
errorbar(dotSize(5), d10.GMADR{i,2}, d10.STDADR{i,2}, 'o', 'color', [0.5 0 0], 'MarkerSize',10, ...
    'MarkerFaceColor', [0.5 0 0])
hold on
errorbar(dotSize(1), d1.GMADR{i,1}, d1.STDADR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(1), d1.GMADL{i,2}, d1.STDADL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
errorbar(dotSize(2), d2.GMADR{i,1}, d2.STDADR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(2), d2.GMADL{i,2}, d2.STDADL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
errorbar(dotSize(3), d3.GMADR{i,1}, d3.STDADR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(3), d3.GMADL{i,2}, d3.STDADL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
errorbar(dotSize(4), d5.GMADR{i,1}, d5.STDADR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(4), d5.GMADL{i,2}, d5.STDADL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
errorbar(dotSize(5), d10.GMADR{i,1}, d10.STDADR{i,1},'o', 'color', 'b','MarkerSize',10, ...
    'MarkerFaceColor','b')
errorbar(dotSize(5), d10.GMADL{i,2}, d10.STDADL{i,2}, 'o', 'color', [0 0 0.5], 'MarkerSize',10, ...
    'MarkerFaceColor', [0 0 0.5])
ylabel('<Vr> FTB - BTF')
xlabel('Dot Size')
set(gca,'xtick',[1 2 3 5 10]);
axis([0 11 0 400])




%%
% ds = [1 2 5 10];
ds = 1;
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);
GMVI = cell(length(ds), 1);
STDVI = cell(length(ds), 1);
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0.5;
    params.AngSize = ds(i);
    %     params.D = 0.1;
    strr = [];
    nn = [];
    vi = [];
    nvi = [];
    for j = 1 : 10
        params.kVisualSize = 0;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -0;
        [outp2] = StartSimulation(params);
        vi = horzcat(vi, (outp2.pss-outp1.pss)/(outp1.pss+outp2.pss));
        nvi = vertcat(nvi, outp1.pss+outp2.pss);
        strr = horzcat(strr, outp1.str);
        nn = vertcat(nn, outp1.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp1.str)])
    end
    GMVI{i} = vi*nvi/sum(nvi);
    STDVI{i} = sqrt((vi-GMVI{i}).*(vi-GMVI{i})*nvi/sum(nvi));
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

% figure,
% errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

figure,
errorbar(ds, cell2mat(GMVI), cell2mat(STDVI)./sqrt(10), 'ok', 'linewidth', 2)
%%

ds = [1 2 5 10];
% ds = 1;
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);
GMVI = cell(length(ds), 1);
STDVI = cell(length(ds), 1);
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0*0.5;
    params.AngSize = ds(i);
    %     params.D = 0.1;
    strr = [];
    nn = [];
    vi = [];
    nvi = [];
    for j = 1 : 10
        params.kVisualSize = 0;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -0;
        [outp2] = StartSimulation(params);
        vi = horzcat(vi, (outp2.pss-outp1.pss)/(outp1.pss+outp2.pss));
        nvi = vertcat(nvi, outp1.pss+outp2.pss);
        strr = horzcat(strr, outp1.str);
        nn = vertcat(nn, outp1.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp1.str)])
    end
    GMVI{i} = vi*nvi/sum(nvi);
    STDVI{i} = sqrt((vi-GMVI{i}).*(vi-GMVI{i})*nvi/sum(nvi));
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

% figure,
% errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

figure,
errorbar(ds, cell2mat(GMVI), cell2mat(STDVI)./sqrt(10), 'ok', 'linewidth', 2)


ds = [1 2 5 10];
% ds = 1;
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);
GMVI = cell(length(ds), 1);
STDVI = cell(length(ds), 1);
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0.5;
    params.AngSize = ds(i);
    %     params.D = 0.1;
    strr = [];
    nn = [];
    vi = [];
    nvi = [];
    for j = 1 : 10
        params.kVisualSize = 80;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -80;
        [outp2] = StartSimulation(params);
        vi = horzcat(vi, (outp2.pss-outp1.pss)/(outp1.pss+outp2.pss));
        nvi = vertcat(nvi, outp1.pss+outp2.pss);
        strr = horzcat(strr, outp1.str);
        nn = vertcat(nn, outp1.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp1.str)])
    end
    GMVI{i} = vi*nvi/sum(nvi);
    STDVI{i} = sqrt((vi-GMVI{i}).*(vi-GMVI{i})*nvi/sum(nvi));
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

% figure,
% errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

figure,
errorbar(ds, cell2mat(GMVI), cell2mat(STDVI)./sqrt(10), 'ok', 'linewidth', 2)


ds = [1 2 5 10];
% ds = 1;
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);
GMVI = cell(length(ds), 1);
STDVI = cell(length(ds), 1);
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0*0.5;
    params.AngSize = ds(i);
    %     params.D = 0.1;
    strr = [];
    nn = [];
    vi = [];
    nvi = [];
    for j = 1 : 10
        params.kVisualSize = 80;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -80;
        [outp2] = StartSimulation(params);
        vi = horzcat(vi, (outp2.pss-outp1.pss)/(outp1.pss+outp2.pss));
        nvi = vertcat(nvi, outp1.pss+outp2.pss);
        strr = horzcat(strr, outp1.str);
        nn = vertcat(nn, outp1.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp1.str)])
    end
    GMVI{i} = vi*nvi/sum(nvi);
    STDVI{i} = sqrt((vi-GMVI{i}).*(vi-GMVI{i})*nvi/sum(nvi));
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

% figure,
% errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

figure,
errorbar(ds, cell2mat(GMVI), cell2mat(STDVI)./sqrt(10), 'ok', 'linewidth', 2)

ds = [1 2 5 10];
% ds = 1;
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);
GMVI = cell(length(ds), 1);
STDVI = cell(length(ds), 1);
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0.5;
    params.AngSize = ds(i);
    %     params.D = 0.1;
    strr = [];
    nn = [];
    vi = [];
    nvi = [];
    for j = 1 : 10
        params.kVisualSize = 90;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -90;
        [outp2] = StartSimulation(params);
        vi = horzcat(vi, (outp2.pss-outp1.pss)/(outp1.pss+outp2.pss));
        nvi = vertcat(nvi, outp1.pss+outp2.pss);
        strr = horzcat(strr, outp1.str);
        nn = vertcat(nn, outp1.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp1.str)])
    end
    GMVI{i} = vi*nvi/sum(nvi);
    STDVI{i} = sqrt((vi-GMVI{i}).*(vi-GMVI{i})*nvi/sum(nvi));
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

% figure,
% errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

figure,
errorbar(ds, cell2mat(GMVI), cell2mat(STDVI)./sqrt(10), 'ok', 'linewidth', 2)


ds = [1 2 5 10];
% ds = 1;
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);
GMVI = cell(length(ds), 1);
STDVI = cell(length(ds), 1);
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0*0.5;
    params.AngSize = ds(i);
    %     params.D = 0.1;
    strr = [];
    nn = [];
    vi = [];
    nvi = [];
    for j = 1 : 10
        params.kVisualSize = 90;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -90;
        [outp2] = StartSimulation(params);
        vi = horzcat(vi, (outp2.pss-outp1.pss)/(outp1.pss+outp2.pss));
        nvi = vertcat(nvi, outp1.pss+outp2.pss);
        strr = horzcat(strr, outp1.str);
        nn = vertcat(nn, outp1.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp1.str)])
    end
    GMVI{i} = vi*nvi/sum(nvi);
    STDVI{i} = sqrt((vi-GMVI{i}).*(vi-GMVI{i})*nvi/sum(nvi));
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

% figure,
% errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

figure,
errorbar(ds, cell2mat(GMVI), cell2mat(STDVI)./sqrt(10), 'ok', 'linewidth', 2)

ds = [1 2 5 10];
% ds = 1;
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);
GMVI = cell(length(ds), 1);
STDVI = cell(length(ds), 1);
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0.5;
    params.AngSize = ds(i);
    %     params.D = 0.1;
    strr = [];
    nn = [];
    vi = [];
    nvi = [];
    for j = 1 : 10
        params.kVisualSize = 100;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -100;
        [outp2] = StartSimulation(params);
        vi = horzcat(vi, (outp2.pss-outp1.pss)/(outp1.pss+outp2.pss));
        nvi = vertcat(nvi, outp1.pss+outp2.pss);
        strr = horzcat(strr, outp1.str);
        nn = vertcat(nn, outp1.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp1.str)])
    end
    GMVI{i} = vi*nvi/sum(nvi);
    STDVI{i} = sqrt((vi-GMVI{i}).*(vi-GMVI{i})*nvi/sum(nvi));
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

% figure,
% errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

figure,
errorbar(ds, cell2mat(GMVI), cell2mat(STDVI)./sqrt(10), 'ok', 'linewidth', 2)


ds = [1 2 5 10];
% ds = 1;
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);
GMVI = cell(length(ds), 1);
STDVI = cell(length(ds), 1);
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0*0.5;
    params.AngSize = ds(i);
    %     params.D = 0.1;
    strr = [];
    nn = [];
    vi = [];
    nvi = [];
    for j = 1 : 10
        params.kVisualSize = 100;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -100;
        [outp2] = StartSimulation(params);
        vi = horzcat(vi, (outp2.pss-outp1.pss)/(outp1.pss+outp2.pss));
        nvi = vertcat(nvi, outp1.pss+outp2.pss);
        strr = horzcat(strr, outp1.str);
        nn = vertcat(nn, outp1.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp1.str)])
    end
    GMVI{i} = vi*nvi/sum(nvi);
    STDVI{i} = sqrt((vi-GMVI{i}).*(vi-GMVI{i})*nvi/sum(nvi));
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

% figure,
% errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

figure,
errorbar(ds, cell2mat(GMVI), cell2mat(STDVI)./sqrt(10), 'ok', 'linewidth', 2)

%% Plots Perfect Agent With Spikes
[params] = GetParams();
params.plotTraj = 1;
params.noisesize = 0;
params.kernelMotorSize = 0;
params.kVisualSize = 0;
StartSimulation(params);



%%
ds = [0.1 1 5 10 15 20 25 30 35];
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = ds(i);
    params.kernelMotorSize = 0*0.5;
    params.kVisualSize = 0*30;
    params.AngSize = 1;
    %     params.D = 0.1;
    strr = [];
    nn = [];
    for j = 1 : 10
        [outp] = StartSimulation(params);
        strr = horzcat(strr, outp.str);
        nn = vertcat(nn, outp.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp.str)])
    end
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

figure,
errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

%%
ds = [0.1 1 10 25 50 75 100];
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0*0.5;
    params.kVisualSize = ds(i);
    params.AngSize = 10;
    %     params.D = 0.1;
    strr = [];
    nn = [];
    for j = 1 : 10
        [outp] = StartSimulation(params);
        strr = horzcat(strr, outp.str);
        nn = vertcat(nn, outp.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp.str)])
    end
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

figure,
errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)

%%
ds = [0.01 0.1 0.5 1 2 5];
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = ds(i);
    params.kVisualSize = 0;
    params.AngSize = 1;
    %     params.D = 0.1;
    strr = [];
    nn = [];
    for j = 1 : 10
        [outp] = StartSimulation(params);
        strr = horzcat(strr, outp.str);
        nn = vertcat(nn, outp.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp.str)])
    end
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

figure,
errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)


%%
ds = [1 2 5 10];
GM = cell(length(ds), 1);
STD = cell(length(ds), 1);

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(length(ds))
end

parfor i = 1:length(ds)
    [params] = GetParams();
    params.L = 6*60*60;
    params.lp_Tau_HR = 15e-3;
    params.hp_Tau_HR = 50e-3;
    params.noisesize = 20;
    params.kernelMotorSize = 0.3;
    params.AngSize = ds(i);
    strr = [];
    nn = [];
    vi = [];
    nvi = [];
    for j = 1 : 10
        params.kVisualSize = 100;
        [outp1] = StartSimulation(params);
        strr = horzcat(strr, outp1.str);
        nn = vertcat(nn, outp1.n);
        disp(['DotSize ' num2str(ds(i)) ', Fly ' num2str(j) ' : ' num2str(outp1.str)])
    end
    GM{i} = strr*nn/sum(nn);
    STD{i} = sqrt((strr-GM{i}).*(strr-GM{i})*nn/sum(nn));
    disp(['DotSize ' num2str(ds(i)) ' : ' num2str(strr*nn/sum(nn)) '+/-' num2str(std(strr))])
end

figure,
errorbar(ds, cell2mat(GM), cell2mat(STD)./sqrt(10), 'ok', 'linewidth', 2)






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% motorWeight = [0 1 2 3 4 5 6 7 8 9 10 12.5 15 17.5 20];
% noiseLevel = [1 2.5 5 7.5 10 15 20 25 30 35 40 45 50 75 100];
motorWeight = 4.9;
noiseLevel = 18;
nSim = 0;
for nL = 1 : length(noiseLevel)
    for mW = 1 : length(motorWeight)
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            parpool(12)
        end
        nC = 12;
        GMM = cell(12,1);
        GN = cell(12,1);
        PSS1M = cell(12,1);
        PSS1N = cell(12,1);
        PSS2M = cell(12,1);
        PSS2N = cell(12,1);
        AL1 = cell(12,1);
        AL2 = cell(12,1);
        PL1 = cell(12,1);
        PL2 = cell(12,1);
        NAL1 = cell(12,1);
        NAL2 = cell(12,1);
        AR1 = cell(12,1);
        AR2 = cell(12,1);
        PR1 = cell(12,1);
        PR2 = cell(12,1);
        NAR1 = cell(12,1);
        NAR2 = cell(12,1);
        parfor k = 1 : 12
            motWei = 4.9;
            noiLev = 18;
            [params] = GetParams();
            params.kernelMotorSize = motWei(mW);
            params.noisesize = noiLev(nL);
            params.kVisualSize = 0;
            params.AngSize = 1;
            params.WeightR = nan;
            params.NCopies = 5;
            [outp1] = StartSimulation(params);
            params.kVisualSize = 0;
            [outp2] = StartSimulation(params);
            GMM{k} = outp1.str;
            GN{k} = outp1.n;
            AL1{k} = outp1.vrl;
            PL1{k} = outp1.pvrl;
            NAL1{k} = outp1.nvrl;
            AL2{k} = outp2.vrl;
            PL2{k} = outp2.pvrl;
            NAL2{k} = outp2.nvrl;
            AR1{k} = outp1.vrr;
            PR1{k} = outp1.pvrr;
            NAR1{k} = outp1.nvrr;
            AR2{k} = outp2.vrr;
            PR2{k} = outp2.pvrr;
            NAR2{k} = outp2.nvrr;
            PSS1M{k} = outp1.pss;
            PSS1N{k} = outp1.npss;
            PSS2M{k} = outp2.pss;
            PSS2N{k} = outp2.npss;
        end
        
        GMM = cell2mat(GMM);
        GN = cell2mat(GN);
        gm = GMM'*GN/sum(GN);
        sem = sqrt(((GMM-gm).*(GMM-gm))'*GN/sum(GN))/sqrt(nC);
        GMS = gm;
        STDS = sem;
        
        GMPSS1M = cell2mat(PSS1M);
        GMPSS1N = cell2mat(PSS1N);
        gmpss1 = GMPSS1M'*GMPSS1N/sum(GMPSS1N);
        sempss1 = sqrt(((GMPSS1M-gmpss1).*(GMPSS1M-gmpss1))'*GMPSS1N/sum(GMPSS1N))/sqrt(nC);
        GMPSS{1} = gmpss1;
        STDPSS{1} = sempss1;
        
        AL1 = cell2mat(AL1);
        NAL1 = cell2mat(NAL1);
        gmal1 = AL1'*NAL1/sum(NAL1);
        semal1 = sqrt(((AL1-gmal1).*(AL1-gmal1))'*NAL1/sum(NAL1))/sqrt(nC);
        GMADL{1} = gmal1;
        STDADL{1} = semal1;
        
        PL1 = cell2mat(PL1);
        gmpl1 = PL1'*NAL1/sum(NAL1);
        sempl1 = sqrt(((PL1-gmpl1).*(PL1-gmpl1))'*NAL1/sum(NAL1))/sqrt(nC);
        GMPL{1} = gmpl1;
        STDPL{1} = sempl1;
        
        AR1 = cell2mat(AR1);
        NAR1 = cell2mat(NAR1);
        gmar1 = AR1'*NAR1/sum(NAR1);
        semar1 = sqrt(((AR1-gmar1).*(AR1-gmar1))'*NAR1/sum(NAR1))/sqrt(nC);
        GMADR{1} = gmar1;
        STDADR{1} = semar1;
        
        PR1 = cell2mat(PR1);
        gmpr1 = PR1'*NAR1/sum(NAR1);
        sempr1 = sqrt(((PR1-gmpr1).*(PR1-gmpr1))'*NAR1/sum(NAR1))/sqrt(nC);
        GMPR{1} = gmpr1;
        STDPR{1} = sempr1;
        
        AL2 = cell2mat(AL2);
        NAL2 = cell2mat(NAL2);
        gmal2 = AL2'*NAL2/sum(NAL2);
        semal2 = sqrt(((AL2-gmal2).*(AL2-gmal2))'*NAL2/sum(NAL2))/sqrt(nC);
        GMADL{2} = gmal2;
        STDADL{2} = semal2;
        
        PL2 = cell2mat(PL2);
        gmpl2 = PL2'*NAL2/sum(NAL2);
        sempl2 = sqrt(((PL2-gmpl2).*(PL2-gmpl2))'*NAL2/sum(NAL2))/sqrt(nC);
        GMPL{2} = gmpl2;
        STDPL{2} = sempl2;
        
        AR2 = cell2mat(AR2);
        NAR2 = cell2mat(NAR2);
        gmar2 = AR2'*NAR2/sum(NAR2);
        semar2 = sqrt(((AR2-gmar2).*(AR2-gmar2))'*NAR2/sum(NAR2))/sqrt(nC);
        GMADR{2} = gmar2;
        STDADR{2} = semar2;
        
        PR2 = cell2mat(PR2);
        gmpr2 = PR2'*NAR2/sum(NAR2);
        sempr2 = sqrt(((PR2-gmpr2).*(PR2-gmpr2))'*NAR2/sum(NAR2))/sqrt(nC);
        GMPR{2} = gmpr2;
        STDPR{2} = sempr2;
        
        GMPSS2M = cell2mat(PSS2M);
        GMPSS2N = cell2mat(PSS2N);
        gmpss2 = GMPSS2M'*GMPSS2N/sum(GMPSS2N);
        sempss2 = sqrt(((GMPSS2M-gmpss2).*(GMPSS2M-gmpss2))'*GMPSS2N/sum(GMPSS2N))/sqrt(10);
        GMPSS{2} = gmpss2;
        STDPSS{2} = sempss2;

        nSim = nSim + 1;
        disp(['Simulation ' num2str(nSim) ': NoiseLevel - ' num2str(nL)  ...
            '; MotorLevel - ' num2str(mW)  ' Str: ' num2str(gm) ...
            ' +/- ' num2str(sem) '  VI:' num2str((gmpss2-gmpss1)/(gmpss2+gmpss1))])
        
        dtt.GMS = GMS;
        dtt.STDS = STDS;
        dtt.GMPSS = GMPSS;
        dtt.STDPSS = STDPSS;
        dtt.GMADL = GMADL;
        dtt.STDADL = STDADL;
        dtt.GMADR = GMADR;
        dtt.STDADR = STDADR;
        dtt.GMPL = GMPL;
        dtt.STDPL = STDPL;
        dtt.GMPR = GMPR;
        dtt.STDPR = STDPR;
        dtt.NL = noiseLevel(nL);
        dtt.MW = motorWeight(mW);
        save(['C:\Users\tomas\Desktop\Simulations Dark\dtDarkNL' num2str(noiseLevel(nL)) ...
            'MW' num2str(motorWeight(mW)) '.mat'], 'dtt');
    end
end



%%

% VisualWeight = [0 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25];
VisualWeight = [0 1 2 3 4 5 6 7 8 9 10 12 14];
% VisualWeight = 6;
% noiseLevel = [10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40];
% noiseLevel = 14.1218;
noiseLevel = [28 30 32 34 36 38 40];

nSim = 0;
for nL = 1 : length(noiseLevel)
    for vW = 1 : length(VisualWeight)
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            parpool(12)
        end
        nC = 12;
        nds = 2;
        GMM = cell(12,nds);
        GN = cell(12,nds);
        PSS1M = cell(12,nds);
        PSS1N = cell(12,nds);
        PSS2M = cell(12,nds);
        PSS2N = cell(12,nds);
        AL1 = cell(12,nds);
        AL2 = cell(12,nds);
        PL1 = cell(12,nds);
        PL2 = cell(12,nds);
        NAL1 = cell(12,nds);
        NAL2 = cell(12,nds);
        AR1 = cell(12,nds);
        AR2 = cell(12,nds);
        PR1 = cell(12,nds);
        PR2 = cell(12,nds);
        NAR1 = cell(12,nds);
        NAR2 = cell(12,nds);
        parfor k = 1 : 12
%             visWei = [0 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25];
                visWei = [0 1 2 3 4 5 6 7 8 9 10 12 14];
%             visWei = 6;
%             noiLev = [10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40];
%             noiLev = 14.1218;
            noiLev = [28 30 32 34 36 38 40];
            DS = [1 2];
            for ds = 1 : nds
                [params] = GetParams();
                params.kernelMotorSize = noiLev(nL)*1.264 - 17.85;
                params.noisesize = noiLev(nL);
                params.kVisualSize = visWei(vW);
                params.AngSize = DS(ds);
                params.WeightR = nan;
                params.NCopies = 3;
                [outp1] = StartSimulation(params);
                params.kVisualSize = -visWei(vW);
                [outp2] = StartSimulation(params);
                GMM{k,ds} = outp1.str;
                GN{k,ds} = outp1.n;
                AL1{k,ds} = outp1.vrl;
                PL1{k,ds} = outp1.pvrl;
                NAL1{k,ds} = outp1.nvrl;
                AL2{k,ds} = outp2.vrl;
                PL2{k,ds} = outp2.pvrl;
                NAL2{k,ds} = outp2.nvrl;
                AR1{k,ds} = outp1.vrr;
                PR1{k,ds} = outp1.pvrr;
                NAR1{k,ds} = outp1.nvrr;
                AR2{k,ds} = outp2.vrr;
                PR2{k,ds} = outp2.pvrr;
                NAR2{k,ds} = outp2.nvrr;
                PSS1M{k,ds} = outp1.pss;
                PSS1N{k,ds} = outp1.npss;
                PSS2M{k,ds} = outp2.pss;
                PSS2N{k,ds} = outp2.npss;
            end
        end
        for ds = 1 : nds           
            if iscell(GMM)
                GMM = cell2mat(GMM);
            end
            if iscell(GN)
                GN = cell2mat(GN);
            end
            gm = GMM(:,ds)'*GN(:,ds)/sum(GN(:,ds));
            sem = sqrt(((GMM(:,ds)-gm).*(GMM(:,ds)-gm))'*...
                GN(:,ds)/sum(GN(:,ds)))/sqrt(nC);
            GMS{ds} = gm;
            STDS{ds} = sem;
            
            if iscell(PSS1M)
                GMPSS1M = cell2mat(PSS1M);
            end
            if iscell(PSS1N)
                GMPSS1N = cell2mat(PSS1N);
            end
            gmpss1 = GMPSS1M(:,ds)'*GMPSS1N(:,ds)/...
                sum(GMPSS1N(:,ds));
            sempss1 = sqrt(((GMPSS1M(:,ds)-gmpss1).*(GMPSS1M(:,ds)-gmpss1))'*...
                GMPSS1N(:,ds)/sum(GMPSS1N(:,ds)))/sqrt(nC);
            GMPSS{1,ds} = gmpss1;
            STDPSS{1,ds} = sempss1;
            
            if iscell(AL1)
                AL1 = cell2mat(AL1);
            end
            if iscell(NAL1)
                NAL1 = cell2mat(NAL1);
            end
            gmal1 = AL1(:,ds)'*NAL1(:,ds)/sum(NAL1(:,ds));
            semal1 = sqrt(((AL1(:,ds)-gmal1).*(AL1(:,ds)-gmal1))'*...
                NAL1(:,ds)/sum(NAL1(:,ds)))/sqrt(nC);
            GMADL{1,ds} = gmal1;
            STDADL{1,ds} = semal1;
            
            if iscell(PL1)
                PL1 = cell2mat(PL1);
            end
            gmpl1 = PL1(:,ds)'*NAL1(:,ds)/sum(NAL1(:,ds));
            sempl1 = sqrt(((PL1(:,ds)-gmpl1).*(PL1(:,ds)-gmpl1))'*...
                NAL1(:,ds)/sum(NAL1(:,ds)))/sqrt(nC);
            GMPL{1,ds} = gmpl1;
            STDPL{1,ds} = sempl1;
            
            if iscell(AR1)
                AR1 = cell2mat(AR1);
            end
            if iscell(NAR1)
                NAR1 = cell2mat(NAR1);
            end
            gmar1 = AR1(:,ds)'*NAR1(:,ds)/sum(NAR1(:,ds));
            semar1 = sqrt(((AR1(:,ds)-gmar1).*(AR1(:,ds)-gmar1))'*...
                NAR1(:,ds)/sum(NAR1(:,ds)))/sqrt(nC);
            GMADR{1,ds} = gmar1;
            STDADR{1,ds} = semar1;
            
            if iscell(PR1)
                PR1 = cell2mat(PR1);
            end
            gmpr1 = PR1(:,ds)'*NAR1(:,ds)/sum(NAR1(:,ds));
            sempr1 = sqrt(((PR1(:,ds)-gmpr1).*(PR1(:,ds)-gmpr1))'*...
                NAR1(:,ds)/sum(NAR1(:,ds)))/sqrt(nC);
            GMPR{1,ds} = gmpr1;
            STDPR{1,ds} = sempr1;
            
            if iscell(AL2)
                AL2 = cell2mat(AL2);
            end
            if iscell(NAL2)
                NAL2 = cell2mat(NAL2);
            end
            gmal2 = AL2(:,ds)'*NAL2(:,ds)/sum(NAL2(:,ds));
            semal2 = sqrt(((AL2(:,ds)-gmal2).*(AL2(:,ds)-gmal2))'*NAL2(:,ds)/sum(NAL2(:,ds)))/sqrt(nC);
            GMADL{2,ds} = gmal2;
            STDADL{2,ds} = semal2;
            
            if iscell(PL2)
                PL2 = cell2mat(PL2);
            end
            gmpl2 = PL2(:,ds)'*NAL2(:,ds)/sum(NAL2(:,ds));
            sempl2 = sqrt(((PL2(:,ds)-gmpl2).*(PL2(:,ds)-gmpl2))'*...
                NAL2(:,ds)/sum(NAL2(:,ds)))/sqrt(nC);
            GMPL{2,ds} = gmpl2;
            STDPL{2,ds} = sempl2;
            
            if iscell(AR2)
                AR2 = cell2mat(AR2);
            end
            if iscell(NAR2)
                NAR2 = cell2mat(NAR2);
            end
            gmar2 = AR2(:,ds)'*NAR2(:,ds)/sum(NAR2(:,ds));
            semar2 = sqrt(((AR2(:,ds)-gmar2).*(AR2(:,ds)-gmar2))'*...
                NAR2(:,ds)/sum(NAR2(:,ds)))/sqrt(nC);
            GMADR{2,ds} = gmar2;
            STDADR{2,ds} = semar2;
            
            if iscell(PR2)
                PR2 = cell2mat(PR2);
            end
            gmpr2 = PR2(:,ds)'*NAR2(:,ds)/sum(NAR2(:,ds));
            sempr2 = sqrt(((PR2(:,ds)-gmpr2).*(PR2(:,ds)-gmpr2))'*...
                NAR2(:,ds)/sum(NAR2(:,ds)))/sqrt(nC);
            GMPR{2,ds} = gmpr2;
            STDPR{2,ds} = sempr2;
            
            if iscell(PSS2M)
                GMPSS2M = cell2mat(PSS2M);
            end
            if iscell(PSS2N)
                GMPSS2N = cell2mat(PSS2N);
            end
            gmpss2 = GMPSS2M(:,ds)'*GMPSS2N(:,ds)/sum(GMPSS2N(:,ds));
            sempss2 = sqrt(((GMPSS2M(:,ds)-gmpss2).*(GMPSS2M(:,ds)-gmpss2))'*...
                GMPSS2N(:,ds)/sum(GMPSS2N(:,ds)))/sqrt(nC);
            GMPSS{2,ds} = gmpss2;
            STDPSS{2,ds} = sempss2;
        end
        
        nSim = nSim + 1;
        disp(['Simulation ' num2str(nSim) ': NoiseLevel - ' num2str(noiseLevel(nL))  ...
            '; MotorLevel - ' num2str(floor(noiseLevel(nL)*1.264 - 17.85)) '; VisualLevel - ' ...
            num2str(VisualWeight(vW))])
        
        dtt.GMS = GMS;
        dtt.STDS = STDS;
        dtt.GMPSS = GMPSS;
        dtt.STDPSS = STDPSS;
        dtt.GMADL = GMADL;
        dtt.STDADL = STDADL;
        dtt.GMADR = GMADR;
        dtt.STDADR = STDADR;
        dtt.GMPL = GMPL;
        dtt.STDPL = STDPL;
        dtt.GMPR = GMPR;
        dtt.STDPR = STDPR;
        dtt.NL = noiseLevel(nL);
        dtt.MW = noiseLevel(nL)*1.264 - 17.85;
        dtt.VW = VisualWeight(vW);
        
        dtt.GMM = GMM;
        dtt.PSS1M = GMPSS1M;
        dtt.PSS2M = GMPSS2M;
        save(['C:\Users\tomas\Desktop\Simulations LightZoom\dtLightNL'...
            num2str(noiseLevel(nL)) 'MW' ...
            num2str(floor(noiseLevel(nL)*1.264 - 17.85)) 'VW'...
            num2str(VisualWeight(vW)) '012_2.mat'], 'dtt');
    end
end
%%
VisualWeight = [0 2 4 6 8 10 15 20 25 30];
MotorWeight = [0 2 4 6 8 10 15 20 25 30 35];
noiseLevel = [15 25 35];

nSim = 0;
for nLL = 1 : length(noiseLevel)
    for mw = 11 : length(MotorWeight)
        for vW = 1 : length(VisualWeight)
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool(12)
            end
            nC = 12;
            nds = 2;
            GMM = cell(12,nds);
            GN = cell(12,nds);
            PSS1M = cell(12,nds);
            PSS1N = cell(12,nds);
            PSS2M = cell(12,nds);
            PSS2N = cell(12,nds);
            AL1 = cell(12,nds);
            AL2 = cell(12,nds);
            PL1 = cell(12,nds);
            PL2 = cell(12,nds);
            NAL1 = cell(12,nds);
            NAL2 = cell(12,nds);
            AR1 = cell(12,nds);
            AR2 = cell(12,nds);
            PR1 = cell(12,nds);
            PR2 = cell(12,nds);
            NAR1 = cell(12,nds);
            NAR2 = cell(12,nds);
            parfor k = 1 : 12
                visW = [0 2 4 6 8 10 15 20 25 30];
                motW = [0 2 4 6 8 10 15 20 25 30 35];
                noiseL = [15 25 35];
                DS = [1 2];
                for ds = 1 : nds
                    [params] = GetParams();
                    params.kernelMotorSize = motW(mw);
                    params.noisesize = noiseL(nLL);
                    params.kVisualSize = visW(vW);
                    params.AngSize = DS(ds);
                    params.WeightR = nan;
                    params.NCopies = 3;
                    [outp1] = StartSimulation(params);
                    params.kVisualSize = -visW(vW);
                    [outp2] = StartSimulation(params);
                    GMM{k,ds} = outp1.str;
                    GN{k,ds} = outp1.n;
                    AL1{k,ds} = outp1.vrl;
                    PL1{k,ds} = outp1.pvrl;
                    NAL1{k,ds} = outp1.nvrl;
                    AL2{k,ds} = outp2.vrl;
                    PL2{k,ds} = outp2.pvrl;
                    NAL2{k,ds} = outp2.nvrl;
                    AR1{k,ds} = outp1.vrr;
                    PR1{k,ds} = outp1.pvrr;
                    NAR1{k,ds} = outp1.nvrr;
                    AR2{k,ds} = outp2.vrr;
                    PR2{k,ds} = outp2.pvrr;
                    NAR2{k,ds} = outp2.nvrr;
                    PSS1M{k,ds} = outp1.pss;
                    PSS1N{k,ds} = outp1.npss;
                    PSS2M{k,ds} = outp2.pss;
                    PSS2N{k,ds} = outp2.npss;
                end
            end
            for ds = 1 : nds
                if iscell(GMM)
                    GMM = cell2mat(GMM);
                end
                if iscell(GN)
                    GN = cell2mat(GN);
                end
                gm = GMM(:,ds)'*GN(:,ds)/sum(GN(:,ds));
                sem = sqrt(((GMM(:,ds)-gm).*(GMM(:,ds)-gm))'*...
                    GN(:,ds)/sum(GN(:,ds)))/sqrt(nC);
                GMS{ds} = gm;
                STDS{ds} = sem;
                
                if iscell(PSS1M)
                    GMPSS1M = cell2mat(PSS1M);
                end
                if iscell(PSS1N)
                    GMPSS1N = cell2mat(PSS1N);
                end
                gmpss1 = GMPSS1M(:,ds)'*GMPSS1N(:,ds)/...
                    sum(GMPSS1N(:,ds));
                sempss1 = sqrt(((GMPSS1M(:,ds)-gmpss1).*(GMPSS1M(:,ds)-gmpss1))'*...
                    GMPSS1N(:,ds)/sum(GMPSS1N(:,ds)))/sqrt(nC);
                GMPSS{1,ds} = gmpss1;
                STDPSS{1,ds} = sempss1;
                
                if iscell(AL1)
                    AL1 = cell2mat(AL1);
                end
                if iscell(NAL1)
                    NAL1 = cell2mat(NAL1);
                end
                gmal1 = AL1(:,ds)'*NAL1(:,ds)/sum(NAL1(:,ds));
                semal1 = sqrt(((AL1(:,ds)-gmal1).*(AL1(:,ds)-gmal1))'*...
                    NAL1(:,ds)/sum(NAL1(:,ds)))/sqrt(nC);
                GMADL{1,ds} = gmal1;
                STDADL{1,ds} = semal1;
                
                if iscell(PL1)
                    PL1 = cell2mat(PL1);
                end
                gmpl1 = PL1(:,ds)'*NAL1(:,ds)/sum(NAL1(:,ds));
                sempl1 = sqrt(((PL1(:,ds)-gmpl1).*(PL1(:,ds)-gmpl1))'*...
                    NAL1(:,ds)/sum(NAL1(:,ds)))/sqrt(nC);
                GMPL{1,ds} = gmpl1;
                STDPL{1,ds} = sempl1;
                
                if iscell(AR1)
                    AR1 = cell2mat(AR1);
                end
                if iscell(NAR1)
                    NAR1 = cell2mat(NAR1);
                end
                gmar1 = AR1(:,ds)'*NAR1(:,ds)/sum(NAR1(:,ds));
                semar1 = sqrt(((AR1(:,ds)-gmar1).*(AR1(:,ds)-gmar1))'*...
                    NAR1(:,ds)/sum(NAR1(:,ds)))/sqrt(nC);
                GMADR{1,ds} = gmar1;
                STDADR{1,ds} = semar1;
                
                if iscell(PR1)
                    PR1 = cell2mat(PR1);
                end
                gmpr1 = PR1(:,ds)'*NAR1(:,ds)/sum(NAR1(:,ds));
                sempr1 = sqrt(((PR1(:,ds)-gmpr1).*(PR1(:,ds)-gmpr1))'*...
                    NAR1(:,ds)/sum(NAR1(:,ds)))/sqrt(nC);
                GMPR{1,ds} = gmpr1;
                STDPR{1,ds} = sempr1;
                
                if iscell(AL2)
                    AL2 = cell2mat(AL2);
                end
                if iscell(NAL2)
                    NAL2 = cell2mat(NAL2);
                end
                gmal2 = AL2(:,ds)'*NAL2(:,ds)/sum(NAL2(:,ds));
                semal2 = sqrt(((AL2(:,ds)-gmal2).*(AL2(:,ds)-gmal2))'*NAL2(:,ds)/sum(NAL2(:,ds)))/sqrt(nC);
                GMADL{2,ds} = gmal2;
                STDADL{2,ds} = semal2;
                
                if iscell(PL2)
                    PL2 = cell2mat(PL2);
                end
                gmpl2 = PL2(:,ds)'*NAL2(:,ds)/sum(NAL2(:,ds));
                sempl2 = sqrt(((PL2(:,ds)-gmpl2).*(PL2(:,ds)-gmpl2))'*...
                    NAL2(:,ds)/sum(NAL2(:,ds)))/sqrt(nC);
                GMPL{2,ds} = gmpl2;
                STDPL{2,ds} = sempl2;
                
                if iscell(AR2)
                    AR2 = cell2mat(AR2);
                end
                if iscell(NAR2)
                    NAR2 = cell2mat(NAR2);
                end
                gmar2 = AR2(:,ds)'*NAR2(:,ds)/sum(NAR2(:,ds));
                semar2 = sqrt(((AR2(:,ds)-gmar2).*(AR2(:,ds)-gmar2))'*...
                    NAR2(:,ds)/sum(NAR2(:,ds)))/sqrt(nC);
                GMADR{2,ds} = gmar2;
                STDADR{2,ds} = semar2;
                
                if iscell(PR2)
                    PR2 = cell2mat(PR2);
                end
                gmpr2 = PR2(:,ds)'*NAR2(:,ds)/sum(NAR2(:,ds));
                sempr2 = sqrt(((PR2(:,ds)-gmpr2).*(PR2(:,ds)-gmpr2))'*...
                    NAR2(:,ds)/sum(NAR2(:,ds)))/sqrt(nC);
                GMPR{2,ds} = gmpr2;
                STDPR{2,ds} = sempr2;
                
                if iscell(PSS2M)
                    GMPSS2M = cell2mat(PSS2M);
                end
                if iscell(PSS2N)
                    GMPSS2N = cell2mat(PSS2N);
                end
                gmpss2 = GMPSS2M(:,ds)'*GMPSS2N(:,ds)/sum(GMPSS2N(:,ds));
                sempss2 = sqrt(((GMPSS2M(:,ds)-gmpss2).*(GMPSS2M(:,ds)-gmpss2))'*...
                    GMPSS2N(:,ds)/sum(GMPSS2N(:,ds)))/sqrt(nC);
                GMPSS{2,ds} = gmpss2;
                STDPSS{2,ds} = sempss2;
            end
            
            nSim = nSim + 1;
            disp(['Simulation ' num2str(nSim) ': NoiseLevel - ' num2str(noiseLevel(nLL))  ...
                '; MotorLevel - ' num2str(MotorWeight(mw)) '; VisualLevel - ' ...
                num2str(VisualWeight(vW))])
            
            dtt.GMS = GMS;
            dtt.STDS = STDS;
            dtt.GMPSS = GMPSS;
            dtt.STDPSS = STDPSS;
            dtt.GMADL = GMADL;
            dtt.STDADL = STDADL;
            dtt.GMADR = GMADR;
            dtt.STDADR = STDADR;
            dtt.GMPL = GMPL;
            dtt.STDPL = STDPL;
            dtt.GMPR = GMPR;
            dtt.STDPR = STDPR;
            dtt.NL = noiseLevel(nLL);
            dtt.MW = MotorWeight(mw);
            dtt.VW = VisualWeight(vW);
            
            dtt.GMM = GMM;
            dtt.PSS1M = GMPSS1M;
            dtt.PSS2M = GMPSS2M;
            save(['C:\Users\tomas\Desktop\SimulationsNL3\dtLightNL'...
                num2str(noiseLevel(nLL)) 'MW' ...
                num2str(MotorWeight(mw)) 'VW'...
                num2str(VisualWeight(vW)) '_12.mat'], 'dtt');
        end
    end
end
