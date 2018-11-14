clear
clc


for dotSize = [1 5 10]
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(12)
    end
    nC = 12;
    nFTBNG = cell(12,1);
    nBTFNG = cell(12,1);
    pssFTBNG = cell(12,1);
    pssBTFNG = cell(12,1);
    vFTBNG = cell(12,1);
    vBTFNG = cell(12,1);
    nFTBRG = cell(12,1);
    nBTFRG = cell(12,1);
    pssFTBRG = cell(12,1);
    pssBTFRG = cell(12,1);
    vFTBRG = cell(12,1);
    vBTFRG = cell(12,1);
    STRNG = cell(12,1);
    NSTRNG = cell(12,1);
    STRRG = cell(12,1);
    NSTRRG = cell(12,1);
    parfor k = 1 : 12
        [params] = GetParams();
        params.kernelMotorSize = 4.9;
        params.noisesize = 15;
        params.kVisualSize = 6.5;
        params.AngSize = dotSize;
        params.WeightR = 1;
        params.NCopies = 5;
        params.L = 10*60*60;
        params.ab = 8.5;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -6.5;
        [outp2] = StartSimulation(params);
        
        STRNG{k} = outp1.str;
        NSTRNG{k} = outp1.n;
        STRRG{k} = outp2.str;
        NSTRRG{k} = outp2.n;
        
        nFTBNG{k} = sum(outp1.nvrr);
        nBTFNG{k} = sum(outp1.nvrl);
        pssFTBNG{k} = outp1.pvrr'*outp1.nvrr/sum(outp1.nvrr);
        pssBTFNG{k} = outp1.pvrl'*outp1.nvrl/sum(outp1.nvrl);
        vFTBNG{k} = outp1.vrr'*outp1.nvrr/sum(outp1.nvrr);
        vBTFNG{k} = outp1.vrl'*outp1.nvrl/sum(outp1.nvrl);
        
        nFTBRG{k} = sum(outp2.nvrl);
        nBTFRG{k} = sum(outp2.nvrr);
        pssFTBRG{k} = outp2.pvrl'*outp2.nvrl/sum(outp2.nvrl);
        pssBTFRG{k} = outp2.pvrr'*outp2.nvrr/sum(outp2.nvrr);
        vFTBRG{k} = outp2.vrl'*outp2.nvrl/sum(outp2.nvrl);
        vBTFRG{k} = outp2.vrr'*outp2.nvrr/sum(outp2.nvrr);
    end
    
    
    NSTRNG = cell2mat(NSTRNG);
    STRNG = cell2mat(STRNG);
    gmSTRNG = STRNG'*NSTRNG/sum(NSTRNG);
    semSTRNG = sqrt(((STRNG-gmSTRNG).*(STRNG-gmSTRNG))'*NSTRNG/sum(NSTRNG))/sqrt(12);
    
    NSTRRG = cell2mat(NSTRRG);
    STRRG = cell2mat(STRRG);
    gmSTRRG = STRRG'*NSTRRG/sum(NSTRNG);
    semSTRRG = sqrt(((STRRG-gmSTRRG).*(STRRG-gmSTRRG))'*NSTRRG/sum(NSTRRG))/sqrt(12);    
    
    
    nFTBNG = cell2mat(nFTBNG);
    gmnFTBNG = mean(nFTBNG);
    semnFTBNG = std(nFTBNG)/sqrt(12);
    
    nFTBRG = cell2mat(nFTBRG);
    gmnFTBRG = mean(nFTBRG);
    semnFTBRG = std(nFTBRG)/sqrt(12);
    
    nBTFNG = cell2mat(nBTFNG);
    gmnBTFNG = mean(nBTFNG);
    semnBTFNG = std(nBTFNG)/sqrt(12);
    
    nBTFRG = cell2mat(nBTFRG);
    gmnBTFRG = mean(nBTFRG);
    semnBTFRG = std(nBTFRG)/sqrt(12);
    
    
    pssFTBNG = cell2mat(pssFTBNG);
    gmpssFTBNG = pssFTBNG'*nFTBNG/sum(nFTBNG);
    sempssFTBNG = sqrt(((pssFTBNG-gmpssFTBNG).*(pssFTBNG-gmpssFTBNG))'*nFTBNG/sum(nFTBNG))/sqrt(12);
    
    pssFTBRG = cell2mat(pssFTBRG);
    gmpssFTBRG = pssFTBRG'*nFTBRG/sum(nFTBRG);
    sempssFTBRG = sqrt(((pssFTBRG-gmpssFTBRG).*(pssFTBRG-gmpssFTBRG))'*nFTBRG/sum(nFTBRG))/sqrt(12);
    
    pssBTFNG = cell2mat(pssBTFNG);
    gmpssBTFNG = pssBTFNG'*nBTFNG/sum(nBTFNG);
    sempssBTFNG = sqrt(((pssBTFNG-gmpssBTFNG).*(pssBTFNG-gmpssBTFNG))'*nBTFNG/sum(nBTFNG))/sqrt(12);
    
    pssBTFRG = cell2mat(pssBTFRG);
    gmpssBTFRG = pssBTFRG'*nBTFRG/sum(nBTFRG);
    sempssBTFRG = sqrt(((pssBTFRG-gmpssBTFRG).*(pssBTFRG-gmpssBTFRG))'*nBTFRG/sum(nBTFRG))/sqrt(12);
    
    
    vFTBNG = cell2mat(vFTBNG);
    gmvFTBNG = vFTBNG'*nFTBNG/sum(nFTBNG);
    semvFTBNG = sqrt(((vFTBNG-gmvFTBNG).*(vFTBNG-gmvFTBNG))'*nFTBNG/sum(nFTBNG))/sqrt(12);
    
    vFTBRG = cell2mat(vFTBRG);
    gmvFTBRG = vFTBRG'*nFTBRG/sum(nFTBRG);
    semvFTBRG = sqrt(((vFTBRG-gmvFTBRG).*(vFTBRG-gmvFTBRG))'*nFTBRG/sum(nFTBRG))/sqrt(12);
    
    vBTFNG = cell2mat(vBTFNG);
    gmvBTFNG = vBTFNG'*nBTFNG/sum(nBTFNG);
    semvBTFNG = sqrt(((vBTFNG-gmvBTFNG).*(vBTFNG-gmvBTFNG))'*nBTFNG/sum(nBTFNG))/sqrt(12);
    
    vBTFRG = cell2mat(vBTFRG);
    gmvBTFRG = vBTFRG'*nBTFRG/sum(nBTFRG);
    semvBTFRG = sqrt(((vBTFRG-gmvBTFRG).*(vBTFRG-gmvBTFRG))'*nBTFRG/sum(nBTFRG))/sqrt(12);
    
    disp('---------------------------------')
    disp('---------------------------------')
    disp(['Dot Size ' num2str(dotSize) ' NG'])
    disp(['nFTB -> ' num2str(gmnFTBNG) '+/-' num2str(semnFTBNG)])
    disp(['nBTF -> ' num2str(gmnBTFNG) '+/-' num2str(semnBTFNG)])
    disp(['pssFTB -> ' num2str(gmpssFTBNG) '+/-' num2str(sempssFTBNG)])
    disp(['pssBTF -> ' num2str(gmpssBTFNG) '+/-' num2str(sempssBTFNG)])
    disp(['vFTB -> ' num2str(gmvFTBNG) '+/-' num2str(semvFTBNG)])
    disp(['vBTF -> ' num2str(gmvBTFNG) '+/-' num2str(semvBTFNG)])
    disp(['str -> ' num2str(gmSTRNG) '+/-' num2str(semSTRNG)])
    disp('---------------------------------')
    disp(['Dot Size ' num2str(dotSize) ' RG'])
    disp(['nFTB -> ' num2str(gmnFTBRG) '+/-' num2str(semnFTBRG)])
    disp(['nBTF -> ' num2str(gmnBTFRG) '+/-' num2str(semnBTFRG)])
    disp(['pssFTB -> ' num2str(gmpssFTBRG) '+/-' num2str(sempssFTBRG)])
    disp(['pssBTF -> ' num2str(gmpssBTFRG) '+/-' num2str(sempssBTFRG)])
    disp(['vFTB -> ' num2str(gmvFTBRG) '+/-' num2str(semvFTBRG)])
    disp(['vBTF -> ' num2str(gmvBTFRG) '+/-' num2str(semvBTFRG)])
    disp(['str -> ' num2str(gmSTRRG) '+/-' num2str(semSTRRG)])
    disp('---------------------------------')
    
    dtt.gmSTRNG = gmSTRNG;
    dtt.semSTRNG = semSTRNG;
    dtt.NSTRNG = NSTRNG;
    dtt.STRNG = STRNG;
    dtt.gmSTRRG = gmSTRRG;
    dtt.semSTRRG = semSTRRG;
    dtt.NSTRRG = NSTRRG;
    dtt.STRRG = STRRG;
    
    dtt.gmnFTBNG = gmnFTBNG;
    dtt.semnFTBNG = semnFTBNG;
    dtt.gmnBTFNG = gmnBTFNG;
    dtt.semnBTFNG = semnBTFNG;
    dtt.gmpssFTBNG = gmpssFTBNG;
    dtt.sempssFTBNG = sempssFTBNG;
    dtt.gmpssBTFNG = gmpssBTFNG;
    dtt.sempssBTFNG = sempssBTFNG;
    dtt.gmvFTBNG = gmvFTBNG;
    dtt.semvFTBNG = semvFTBNG;
    dtt.gmvBTFNG = gmvBTFNG;
    dtt.semvBTFNG = semvBTFNG;
    dtt.nFTBNG = nFTBNG;
    dtt.nBTFNG = nBTFNG;
    dtt.pssFTBNG = pssFTBNG;
    dtt.pssBTFNG = pssBTFNG;
    dtt.vFTBNG = vFTBNG;
    dtt.vBTFNG = vBTFNG;
    
    dtt.gmnFTBRG = gmnFTBRG;
    dtt.semnFTBRG = semnFTBRG;
    dtt.gmnBTFRG = gmnBTFRG;
    dtt.semnBTFRG = semnBTFRG;
    dtt.gmpssFTBRG = gmpssFTBRG;
    dtt.sempssFTBRG = sempssFTBRG;
    dtt.gmpssBTFRG = gmpssBTFRG;
    dtt.sempssBTFRG = sempssBTFRG;
    dtt.gmvFTBRG = gmvFTBRG;
    dtt.semvFTBRG = semvFTBRG;
    dtt.gmvBTFRG = gmvBTFRG;
    dtt.semvBTFRG = semvBTFRG;
    dtt.nFTBRG = nFTBRG;
    dtt.nBTFRG = nBTFRG;
    dtt.pssFTBRG = pssFTBRG;
    dtt.pssBTFRG = pssBTFRG;
    dtt.vFTBRG = vFTBRG;
    dtt.vBTFRG = vBTFRG;
    
    save(['C:\Users\tomas\Desktop\SimAB8.5' num2str(dotSize) '.mat'], 'dtt');
end



clear
clc
%%

for dotSize = [1 5 10]
%     dotSize = 10;
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(12)
    end
    nC = 12;
    nFTBNG = cell(12,1);
    nBTFNG = cell(12,1);
    pssFTBNG = cell(12,1);
    pssBTFNG = cell(12,1);
    vFTBNG = cell(12,1);
    vBTFNG = cell(12,1);
    nFTBRG = cell(12,1);
    nBTFRG = cell(12,1);
    pssFTBRG = cell(12,1);
    pssBTFRG = cell(12,1);
    vFTBRG = cell(12,1);
    vBTFRG = cell(12,1);
    STRNG = cell(12,1);
    NSTRNG = cell(12,1);
    STRRG = cell(12,1);
    NSTRRG = cell(12,1);
    parfor k = 1 : 12
        [params] = GetParams();
        params.kernelMotorSize = 4.9;
        params.noisesize = 15;
        params.kVisualSize = 6.5;
        params.AngSize = dotSize;
        params.WeightR = 1;
        params.NCopies = 5;
        params.L = 10*60*60;
%         params.ab = 12.36;
        [outp1] = StartSimulation(params);
        params.kVisualSize = -6.5;
        [outp2] = StartSimulation(params);
        
        STRNG{k} = outp1.str;
        NSTRNG{k} = outp1.n;
        STRRG{k} = outp2.str;
        NSTRRG{k} = outp2.n;
        
        nFTBNG{k} = sum(outp1.nvrr);
        nBTFNG{k} = sum(outp1.nvrl);
        pssFTBNG{k} = outp1.pvrr'*outp1.nvrr/sum(outp1.nvrr);
        pssBTFNG{k} = outp1.pvrl'*outp1.nvrl/sum(outp1.nvrl);
        vFTBNG{k} = outp1.vrr'*outp1.nvrr/sum(outp1.nvrr);
        vBTFNG{k} = outp1.vrl'*outp1.nvrl/sum(outp1.nvrl);
        
        nFTBRG{k} = sum(outp2.nvrl);
        nBTFRG{k} = sum(outp2.nvrr);
        pssFTBRG{k} = outp2.pvrl'*outp2.nvrl/sum(outp2.nvrl);
        pssBTFRG{k} = outp2.pvrr'*outp2.nvrr/sum(outp2.nvrr);
        vFTBRG{k} = outp2.vrl'*outp2.nvrl/sum(outp2.nvrl);
        vBTFRG{k} = outp2.vrr'*outp2.nvrr/sum(outp2.nvrr);
    end
    
    NSTRNG = cell2mat(NSTRNG);
    STRNG = cell2mat(STRNG);
    gmSTRNG = STRNG'*NSTRNG/sum(NSTRNG);
    semSTRNG = sqrt(((STRNG-gmSTRNG).*(STRNG-gmSTRNG))'*NSTRNG/sum(NSTRNG))/sqrt(12);
    
    NSTRRG = cell2mat(NSTRRG);
    STRRG = cell2mat(STRRG);
    gmSTRRG = STRRG'*NSTRRG/sum(NSTRNG);
    semSTRRG = sqrt(((STRRG-gmSTRRG).*(STRRG-gmSTRRG))'*NSTRRG/sum(NSTRRG))/sqrt(12);    
    
    
    nFTBNG = cell2mat(nFTBNG);
    gmnFTBNG = mean(nFTBNG);
    semnFTBNG = std(nFTBNG)/sqrt(12);
    
    nFTBRG = cell2mat(nFTBRG);
    gmnFTBRG = mean(nFTBRG);
    semnFTBRG = std(nFTBRG)/sqrt(12);
    
    nBTFNG = cell2mat(nBTFNG);
    gmnBTFNG = mean(nBTFNG);
    semnBTFNG = std(nBTFNG)/sqrt(12);
    
    nBTFRG = cell2mat(nBTFRG);
    gmnBTFRG = mean(nBTFRG);
    semnBTFRG = std(nBTFRG)/sqrt(12);
    
    
    pssFTBNG = cell2mat(pssFTBNG);
    gmpssFTBNG = pssFTBNG'*nFTBNG/sum(nFTBNG);
    sempssFTBNG = sqrt(((pssFTBNG-gmpssFTBNG).*(pssFTBNG-gmpssFTBNG))'*nFTBNG/sum(nFTBNG))/sqrt(12);
    
    pssFTBRG = cell2mat(pssFTBRG);
    gmpssFTBRG = pssFTBRG'*nFTBRG/sum(nFTBRG);
    sempssFTBRG = sqrt(((pssFTBRG-gmpssFTBRG).*(pssFTBRG-gmpssFTBRG))'*nFTBRG/sum(nFTBRG))/sqrt(12);
    
    pssBTFNG = cell2mat(pssBTFNG);
    gmpssBTFNG = pssBTFNG'*nBTFNG/sum(nBTFNG);
    sempssBTFNG = sqrt(((pssBTFNG-gmpssBTFNG).*(pssBTFNG-gmpssBTFNG))'*nBTFNG/sum(nBTFNG))/sqrt(12);
    
    pssBTFRG = cell2mat(pssBTFRG);
    gmpssBTFRG = pssBTFRG'*nBTFRG/sum(nBTFRG);
    sempssBTFRG = sqrt(((pssBTFRG-gmpssBTFRG).*(pssBTFRG-gmpssBTFRG))'*nBTFRG/sum(nBTFRG))/sqrt(12);
    
    
    vFTBNG = cell2mat(vFTBNG);
    gmvFTBNG = vFTBNG'*nFTBNG/sum(nFTBNG);
    semvFTBNG = sqrt(((vFTBNG-gmvFTBNG).*(vFTBNG-gmvFTBNG))'*nFTBNG/sum(nFTBNG))/sqrt(12);
    
    vFTBRG = cell2mat(vFTBRG);
    gmvFTBRG = vFTBRG'*nFTBRG/sum(nFTBRG);
    semvFTBRG = sqrt(((vFTBRG-gmvFTBRG).*(vFTBRG-gmvFTBRG))'*nFTBRG/sum(nFTBRG))/sqrt(12);
    
    vBTFNG = cell2mat(vBTFNG);
    gmvBTFNG = vBTFNG'*nBTFNG/sum(nBTFNG);
    semvBTFNG = sqrt(((vBTFNG-gmvBTFNG).*(vBTFNG-gmvBTFNG))'*nBTFNG/sum(nBTFNG))/sqrt(12);
    
    vBTFRG = cell2mat(vBTFRG);
    gmvBTFRG = vBTFRG'*nBTFRG/sum(nBTFRG);
    semvBTFRG = sqrt(((vBTFRG-gmvBTFRG).*(vBTFRG-gmvBTFRG))'*nBTFRG/sum(nBTFRG))/sqrt(12);
    
    disp('---------------------------------')
    disp('---------------------------------')
    disp(['Dot Size ' num2str(dotSize) ' NG'])
    disp(['nFTB -> ' num2str(gmnFTBNG) '+/-' num2str(semnFTBNG)])
    disp(['nBTF -> ' num2str(gmnBTFNG) '+/-' num2str(semnBTFNG)])
    disp(['pssFTB -> ' num2str(gmpssFTBNG) '+/-' num2str(sempssFTBNG)])
    disp(['pssBTF -> ' num2str(gmpssBTFNG) '+/-' num2str(sempssBTFNG)])
    disp(['vFTB -> ' num2str(gmvFTBNG) '+/-' num2str(semvFTBNG)])
    disp(['vBTF -> ' num2str(gmvBTFNG) '+/-' num2str(semvBTFNG)])
    disp(['str -> ' num2str(gmSTRNG) '+/-' num2str(semSTRNG)])
    disp('---------------------------------')
    disp(['Dot Size ' num2str(dotSize) ' RG'])
    disp(['nFTB -> ' num2str(gmnFTBRG) '+/-' num2str(semnFTBRG)])
    disp(['nBTF -> ' num2str(gmnBTFRG) '+/-' num2str(semnBTFRG)])
    disp(['pssFTB -> ' num2str(gmpssFTBRG) '+/-' num2str(sempssFTBRG)])
    disp(['pssBTF -> ' num2str(gmpssBTFRG) '+/-' num2str(sempssBTFRG)])
    disp(['vFTB -> ' num2str(gmvFTBRG) '+/-' num2str(semvFTBRG)])
    disp(['vBTF -> ' num2str(gmvBTFRG) '+/-' num2str(semvBTFRG)])
    disp(['str -> ' num2str(gmSTRRG) '+/-' num2str(semSTRRG)])
    disp('---------------------------------')
    
    dtt.gmSTRNG = gmSTRNG;
    dtt.semSTRNG = semSTRNG;
    dtt.NSTRNG = NSTRNG;
    dtt.STRNG = STRNG;
    dtt.gmSTRRG = gmSTRRG;
    dtt.semSTRRG = semSTRRG;
    dtt.NSTRRG = NSTRRG;
    dtt.STRRG = STRRG;
    
    
    dtt.gmnFTBNG = gmnFTBNG;
    dtt.semnFTBNG = semnFTBNG;
    dtt.gmnBTFNG = gmnBTFNG;
    dtt.semnBTFNG = semnBTFNG;
    dtt.gmpssFTBNG = gmpssFTBNG;
    dtt.sempssFTBNG = sempssFTBNG;
    dtt.gmpssBTFNG = gmpssBTFNG;
    dtt.sempssBTFNG = sempssBTFNG;
    dtt.gmvFTBNG = gmvFTBNG;
    dtt.semvFTBNG = semvFTBNG;
    dtt.gmvBTFNG = gmvBTFNG;
    dtt.semvBTFNG = semvBTFNG;
    dtt.nFTBNG = nFTBNG;
    dtt.nBTFNG = nBTFNG;
    dtt.pssFTBNG = pssFTBNG;
    dtt.pssBTFNG = pssBTFNG;
    dtt.vFTBNG = vFTBNG;
    dtt.vBTFNG = vBTFNG;
    
    dtt.gmnFTBRG = gmnFTBRG;
    dtt.semnFTBRG = semnFTBRG;
    dtt.gmnBTFRG = gmnBTFRG;
    dtt.semnBTFRG = semnBTFRG;
    dtt.gmpssFTBRG = gmpssFTBRG;
    dtt.sempssFTBRG = sempssFTBRG;
    dtt.gmpssBTFRG = gmpssBTFRG;
    dtt.sempssBTFRG = sempssBTFRG;
    dtt.gmvFTBRG = gmvFTBRG;
    dtt.semvFTBRG = semvFTBRG;
    dtt.gmvBTFRG = gmvBTFRG;
    dtt.semvBTFRG = semvBTFRG;
    dtt.nFTBRG = nFTBRG;
    dtt.nBTFRG = nBTFRG;
    dtt.pssFTBRG = pssFTBRG;
    dtt.pssBTFRG = pssBTFRG;
    dtt.vFTBRG = vFTBRG;
    dtt.vBTFRG = vBTFRG;
    
%     save(['C:\Users\tomas\Desktop\SimAB12.36_2' num2str(dotSize) '.mat'], 'dtt');
end