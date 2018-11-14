clear
clc
path = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Software\Simulations\Simulations3NL\';
cnds = dir(path);
cnds = cnds(3:end);
str = cell(5,1);
errstr = cell(5,1);
vi = cell(5,1);
errvi = cell(5,1);
nl = cell(5,1);
mw = cell(5,1);
vw = cell(5,1);
for i = 1 : length(cnds)
    dt = load([path cnds(i).name]);
    dt = dt.dtt;
    if length(dt.GMS)==3
        str{3} = vertcat(str{3}, dt.GMS{1});
        str{4} = vertcat(str{4}, dt.GMS{2});
        str{5} = vertcat(str{5}, dt.GMS{3});
        errstr{3} = vertcat(errstr{3}, dt.STDS{1});
        errstr{4} = vertcat(errstr{4}, dt.STDS{2});
        errstr{5} = vertcat(errstr{5}, dt.STDS{3});
        viaux = (dt.GMPSS{2,1}-dt.GMPSS{1,1})/(dt.GMPSS{1,1}+dt.GMPSS{2,1});
        errviaux = sqrt((dt.STDPSS{2,1}*2*dt.GMPSS{1,1}/((dt.GMPSS{1,1}+dt.GMPSS{2,1})))^2 ...
            + (dt.STDPSS{1,1}*2*dt.GMPSS{2,1}/((dt.GMPSS{1,1}+dt.GMPSS{2,1})))^2);
        vi{3} = vertcat(vi{3}, viaux);
        errvi{3} = vertcat(errvi{3}, errviaux);
        viaux = (dt.GMPSS{2,2}-dt.GMPSS{1,2})/(dt.GMPSS{1,2}+dt.GMPSS{2,2});
        errviaux = sqrt((dt.STDPSS{2,2}*2*dt.GMPSS{1,2}/((dt.GMPSS{1,2}+dt.GMPSS{2,2})))^2 ...
            + (dt.STDPSS{1,2}*2*dt.GMPSS{2,2}/((dt.GMPSS{1,2}+dt.GMPSS{2,2})))^2);
        vi{4} = vertcat(vi{4}, viaux);
        errvi{4} = vertcat(errvi{4}, errviaux);
        viaux = (dt.GMPSS{2,3}-dt.GMPSS{1,3})/(dt.GMPSS{1,3}+dt.GMPSS{2,3});
        errviaux = sqrt((dt.STDPSS{2,3}*2*dt.GMPSS{1,3}/((dt.GMPSS{1,3}+dt.GMPSS{2,3})))^2 ...
            + (dt.STDPSS{1,3}*2*dt.GMPSS{2,3}/((dt.GMPSS{1,3}+dt.GMPSS{2,3})))^2);
        vi{5} = vertcat(vi{5}, viaux);
        errvi{5} = vertcat(errvi{5}, errviaux);
        nl{3} = vertcat(nl{3}, dt.NL);
        mw{3} = vertcat(mw{3}, dt.MW);
        vw{3} = vertcat(vw{3}, dt.VW);
        nl{4} = vertcat(nl{4}, dt.NL);
        mw{4} = vertcat(mw{4}, dt.MW);
        vw{4} = vertcat(vw{4}, dt.VW);
        nl{5} = vertcat(nl{5}, dt.NL);
        mw{5} = vertcat(mw{5}, dt.MW);
        vw{5} = vertcat(vw{5}, dt.VW);
    else
        str{1} = vertcat(str{1}, dt.GMS{1});
        str{2} = vertcat(str{2}, dt.GMS{2});
        errstr{1} = vertcat(errstr{1}, dt.STDS{1});
        errstr{2} = vertcat(errstr{2}, dt.STDS{2});
        viaux = (dt.GMPSS{2,1}-dt.GMPSS{1,1})/(dt.GMPSS{1,1}+dt.GMPSS{2,1});
        errviaux = sqrt((dt.STDPSS{2,1}*2*dt.GMPSS{1,1}/((dt.GMPSS{1,1}+dt.GMPSS{2,1})))^2 ...
            + (dt.STDPSS{1,1}*2*dt.GMPSS{2,1}/((dt.GMPSS{1,1}+dt.GMPSS{2,1})))^2);
        vi{1} = vertcat(vi{1}, viaux);
        errvi{1} = vertcat(errvi{1}, errviaux);
        viaux = (dt.GMPSS{2,2}-dt.GMPSS{1,2})/(dt.GMPSS{1,2}+dt.GMPSS{2,2});
        errviaux = sqrt((dt.STDPSS{2,2}*2*dt.GMPSS{1,2}/((dt.GMPSS{1,2}+dt.GMPSS{2,2})))^2 ...
            + (dt.STDPSS{1,2}*2*dt.GMPSS{2,2}/((dt.GMPSS{1,2}+dt.GMPSS{2,2})))^2);
        vi{2} = vertcat(vi{2}, viaux);
        errvi{2} = vertcat(errvi{2}, errviaux);
        nl{1} = vertcat(nl{1}, dt.NL);
        mw{1} = vertcat(mw{1}, dt.MW);
        vw{1} = vertcat(vw{1}, dt.VW);
        nl{2} = vertcat(nl{2}, dt.NL);
        mw{2} = vertcat(mw{2}, dt.MW);
        vw{2} = vertcat(vw{2}, dt.VW);
    end
end

avdt(1) = 31;
avdt(2) = 36.5;
avdt(3) = 41;
avdt(4) = 45;
avdt(5) = 47.5;
avdt(6) = 52.5;
AvErr = 1.5;


avdtr(1) = 31;
avdtr(2) = 36.5;
avdtr(3) = 43;
avdtr(4) = 47.5;
avdtr(5) = 52.5;
AvErr = 1.5;

%%
VW = [0 2 4 6 8 10 15 20 25 30];
MW = [0 2 4 6 8 10 15 20 25 30];
NL = [15 25 35];
ERR = cell(3,1);
STR = cell(3,1);
ESTR = cell(3,1+);
for nnl = 1 : 3
    ERR{nnl} = zeros(10,10,1);
    STR{nnl} = zeros(10,10,6);
    ESTR{nnl} = zeros(10,10,6);
    for mmw = 1 : 10
        ind = find(vw{1} == 0 & nl{1} == NL(nnl) & mw{1} == MW(mmw));
        err0 = 0;
        str0 = 0;
        estr0 = 0;
        for i = 1 : 5
            err0 = err0 + (str{i}(ind)-avdt(1))^2;
            str0 = str0 + str{i}(ind);
            estr0 = estr0 + errstr{i}(ind);
        end
        str0 = str0/5;
        estr0 = estr0/5;
        for vvw = 1 : 10
            ind = find(vw{1} == VW(vvw) & nl{1} == NL(nnl) & mw{1} == MW(mmw));
            STR{nnl}(mmw, vvw, 1) = str0;
            ESTR{nnl}(mmw, vvw, 1) = estr0;
            errST = 0;
            for i = 1 : 5
                errST = errST + (str{i}(ind)-avdt(1+i))^2;
                STR{nnl}(mmw, vvw, 1+i) = str{i}(ind);
                ESTR{nnl}(mmw, vvw, 1+i) = errstr{i}(ind);
            end
            errST = sqrt(errST +err0);
            ERR{nnl}(mmw, vvw, 1) = errST;
        end
    end
end

%%
minM = min(ERR{1}(:));
[rowM(1),colM(1)] = find(ERR{1}==minM);
minM = min(ERR{2}(:));
[rowM(2),colM(2)] = find(ERR{2}==minM);
minM = min(ERR{3}(:));
[rowM(3),colM(3)] = find(ERR{3}==minM);

[MinNM(1), MinNMI(1)] = min(ERR{1}(1,:));
[MinNM(2), MinNMI(2)] = min(ERR{2}(1,:));
[MinNM(3), MinNMI(3)] = min(ERR{3}(1,:));

[MinNV(1), MinNVI(1)] = min(ERR{1}(:,1));
[MinNV(2), MinNVI(2)] = min(ERR{2}(:,1));
[MinNV(3), MinNVI(3)] = min(ERR{3}(:,1));



NLS = [15 25 35];
ds = [0 1 2 3 5 10];
ds2 = [1 2.5 5 10];
cmap = autumn(3);
cmap3 = summer(4);
figure,
hold on
for n = 1 : 3
    subplot(2,3,n)
    hold on
    errorbar(0, avdtr(1), AvErr, 'ok')
    for p = 1 : 4
%         errorbar(ds2(p), avdtr(p+1), AvErr, 'o', 'color', cmap2(p,:))
        errorbar(ds2(p), avdtr(p+1), AvErr, 'o', 'color', 'b')
    end
    hold on
%     for i = 1 : 6
%         errorbar(ds(i), STR{n}(rowM(n),colM(n),i), ESTR{n}(rowM(n),colM(n),i), 'o', 'color', cmap(n,:))
%         axis([-1 11 10 70])
%     end
%     for i = 1 : 6
%         errorbar(ds(i), STR{n}(1,MinNMI(n),i), ESTR{n}(1,MinNMI(n),i), 'o', 'color', cmap3(4-n,:))
%         axis([-1 11 10 70])
%     end
    for i = 1 : 6
        errorbar(ds(i), STR{n}(MinNVI(n),1,i), ESTR{n}(MinNVI(n),1,i), 'o', 'color', cmap3(4-n,:))
        axis([-1 11 10 70])
    end
    title(['NL = ' num2str(NLS(n)) ', vw = ' num2str(VW(colM(n))) ...
        ', mw = ' num2str(MW(rowM(n))) ])
    xlabel('Dot Size (º)')
    ylabel('Straightness (a.u.)')
    
    subplot(2,3,[4 5 6])
    bar(n, ERR{n}(rowM(n),colM(n)),'FaceColor', cmap(n,:))
    hold on
    bar(4+n, ERR{n}(1,MinNMI(n)),'FaceColor', cmap(n,:))
    hold on
    bar(8+n, ERR{n}(MinNVI(n),1), 'FaceColor',cmap(n,:))
    xlabel('Condition')
    ylabel('Min Error (a.u.)')
    xt={'V+M' ; 'V' ; 'M'} ;
    set(gca,'xtick',[2 6 10]);
    set(gca,'xticklabel',xt);
end












%% PerFLy
clear
clc
path = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Software\Simulations\Simulations3NL\';
cnds = dir(path);
cnds = cnds(3:end);
str = cell(5,12);
errstr = cell(5,12);
vi = cell(5,12);
errvi = cell(5,12);
nl = cell(5,12);
mw = cell(5,12);
vw = cell(5,12);
for i = 1 : length(cnds)
    dt = load([path cnds(i).name]);
    dt = dt.dtt;
    if length(dt.GMS) == 3

        
        str{3} = vertcat(str{3}, dt.GMS{1});
        str{4} = vertcat(str{4}, dt.GMS{2});
        str{5} = vertcat(str{5}, dt.GMS{3});
        errstr{3} = vertcat(errstr{3}, dt.STDS{1});
        errstr{4} = vertcat(errstr{4}, dt.STDS{2});
        errstr{5} = vertcat(errstr{5}, dt.STDS{3});
        viaux = (dt.GMPSS{2,1}-dt.GMPSS{1,1})/(dt.GMPSS{1,1}+dt.GMPSS{2,1});
        errviaux = sqrt((dt.STDPSS{2,1}*2*dt.GMPSS{1,1}/((dt.GMPSS{1,1}+dt.GMPSS{2,1})))^2 ...
            + (dt.STDPSS{1,1}*2*dt.GMPSS{2,1}/((dt.GMPSS{1,1}+dt.GMPSS{2,1})))^2);
        vi{3} = vertcat(vi{3}, viaux);
        errvi{3} = vertcat(errvi{3}, errviaux);
        viaux = (dt.GMPSS{2,2}-dt.GMPSS{1,2})/(dt.GMPSS{1,2}+dt.GMPSS{2,2});
        errviaux = sqrt((dt.STDPSS{2,2}*2*dt.GMPSS{1,2}/((dt.GMPSS{1,2}+dt.GMPSS{2,2})))^2 ...
            + (dt.STDPSS{1,2}*2*dt.GMPSS{2,2}/((dt.GMPSS{1,2}+dt.GMPSS{2,2})))^2);
        vi{4} = vertcat(vi{4}, viaux);
        errvi{4} = vertcat(errvi{4}, errviaux);
        viaux = (dt.GMPSS{2,3}-dt.GMPSS{1,3})/(dt.GMPSS{1,3}+dt.GMPSS{2,3});
        errviaux = sqrt((dt.STDPSS{2,3}*2*dt.GMPSS{1,3}/((dt.GMPSS{1,3}+dt.GMPSS{2,3})))^2 ...
            + (dt.STDPSS{1,3}*2*dt.GMPSS{2,3}/((dt.GMPSS{1,3}+dt.GMPSS{2,3})))^2);
        vi{5} = vertcat(vi{5}, viaux);
        errvi{5} = vertcat(errvi{5}, errviaux);
        nl{3} = vertcat(nl{3}, dt.NL);
        mw{3} = vertcat(mw{3}, dt.MW);
        vw{3} = vertcat(vw{3}, dt.VW);
        nl{4} = vertcat(nl{4}, dt.NL);
        mw{4} = vertcat(mw{4}, dt.MW);
        vw{4} = vertcat(vw{4}, dt.VW);
        nl{5} = vertcat(nl{5}, dt.NL);
        mw{5} = vertcat(mw{5}, dt.MW);
        vw{5} = vertcat(vw{5}, dt.VW);
    else
        str{1} = vertcat(str{1}, dt.GMS{1});
        str{2} = vertcat(str{2}, dt.GMS{2});
        errstr{1} = vertcat(errstr{1}, dt.STDS{1});
        errstr{2} = vertcat(errstr{2}, dt.STDS{2});
        viaux = (dt.GMPSS{2,1}-dt.GMPSS{1,1})/(dt.GMPSS{1,1}+dt.GMPSS{2,1});
        errviaux = sqrt((dt.STDPSS{2,1}*2*dt.GMPSS{1,1}/((dt.GMPSS{1,1}+dt.GMPSS{2,1})))^2 ...
            + (dt.STDPSS{1,1}*2*dt.GMPSS{2,1}/((dt.GMPSS{1,1}+dt.GMPSS{2,1})))^2);
        vi{1} = vertcat(vi{1}, viaux);
        errvi{1} = vertcat(errvi{1}, errviaux);
        viaux = (dt.GMPSS{2,2}-dt.GMPSS{1,2})/(dt.GMPSS{1,2}+dt.GMPSS{2,2});
        errviaux = sqrt((dt.STDPSS{2,2}*2*dt.GMPSS{1,2}/((dt.GMPSS{1,2}+dt.GMPSS{2,2})))^2 ...
            + (dt.STDPSS{1,2}*2*dt.GMPSS{2,2}/((dt.GMPSS{1,2}+dt.GMPSS{2,2})))^2);
        vi{2} = vertcat(vi{2}, viaux);
        errvi{2} = vertcat(errvi{2}, errviaux);
        nl{1} = vertcat(nl{1}, dt.NL);
        mw{1} = vertcat(mw{1}, dt.MW);
        vw{1} = vertcat(vw{1}, dt.VW);
        nl{2} = vertcat(nl{2}, dt.NL);
        mw{2} = vertcat(mw{2}, dt.MW);
        vw{2} = vertcat(vw{2}, dt.VW);
    end
end



