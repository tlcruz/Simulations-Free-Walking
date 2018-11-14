clear
clc
% path = 'C:\Users\tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Software\Simulations\\SimulationsLight\';
path = 'C:\Users\Tomas\Dropbox (Sensorimotor)\ChiappeLabNew\data\TOMÁS\Free Walking VR\Software\Simulations\Simulations No Motor\';
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


%%
avdt(1) = 36.5;
avdt(2) = 41;
avdt(3) = 45;
avdt(4) = 47.5;
avdt(5) = 52.5;

ds = 0;%[1 2 3 5 10];

minE = 5;
cmap = jet(minE);
X = [];
Y = [];
Z = [];
figure,
hold on
errVect = [];
i = 1;
for j = 1 : length(str{i})
    errST = 0;
    for i = 1 : 5
        errST = errST + (str{i}(j)-avdt(i))^2;
    end
    errST = sqrt(errST);
    errVect = vertcat(errVect,errST);
    X = vertcat(X, vw{i}(j));
    Y = vertcat(Y, nl{i}(j));
    Z = vertcat(Z, errST);
    scatter(vw{i}(j),nl{i}(j), 100,  cmap(min(minE,floor(errST)),:), 'filled')
end
axis([0 25 12 34])
%%
figure,
[xq,yq] = meshgrid(0:0.1:100, 10:0.1:35);
zq = griddata(X,Y,Z,xq,yq);
s = surf(xq,yq,zq,'FaceAlpha',0.7);
% surf(xq,yq,zq)
% hold on
% plot3(X,Y,Z,'o')
colormap jet
s.EdgeColor = 'none';
xlim([0 25])
ylim([10 35])
% axis off
% xlabel('Visual Weight (a.u.)')
% ylabel('Noise Level (a.u.)')
% zlabel('NDDist (a.u.)')

%%
i = 1;
% noiseLevel = [10 12 14 16 18 20 22 24 26 28 30 32 34 36 38];
noiseLevel = [14 16 18 20 22 24 26 28 30 32];
figure,
hold on
for n = 1 : length(noiseLevel)
    auxErr = [];
    for j = 1 : length(Y)
        if Y(j) == noiseLevel(n)
            auxErr = vertcat(auxErr, Z(j));
        end
    end
    
    scatter(noiseLevel(n), min(auxErr), 100, 'k', 'filled')
    axis([9 35 0 10])
end
%%
[m,I] = sort(errVect);
I = I(2);
figure,
cmap = autumn(5);
ds = [1 2 3 5 10];
subplot(1,2,1)

for i = 1 : 5
%     subplot(1,5,i)
    hold on
    for j = 1 : length(str{i})
        if nl{i}(j) == Y(I) && vw{i}(j) == X(I)
            errorbar(ds(i), vi{i}(j), errvi{i}(j),...
                'o', 'color', cmap(i,:))
            axis([0 11 -0.05 0.35])
        end
    end
    title(['Dot Size: ' num2str(ds(i))])
end

% figure, 
subplot(1,2,2)

hold on
for i = 1 : 5
    for j = 1 : length(str{i})
        if nl{i}(j) == Y(I) && vw{i}(j) == X(I)
            errorbar(ds(i), str{i}(j), errstr{i}(j), 'o', 'color', cmap(i,:))
            axis([0 11 20 60])
        end
    end
end

%%
figure,
hold on
for i = 1 : 5
    for j = 1 : length(str{i})
        if nl{i}(j) == Y(I) && vw{i}(j) == X(I)
            errorbar(ds(i), str{i}(j), str{i}(j), cmap(i,:), 'filled')
            axis([0 11 20 60])
        end
    end
end

%%
noiseLevel = [10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40];
for ds = 1 : 5
    for i = 1 : length(noiseLevel)
        disp([num2str(ds) '->   ' num2str(noiseLevel(i)) ': ' num2str(length(find(nl{ds}==noiseLevel(i))==1))])
    end
end
disp('   ')
























%%
maxdt(1) = 38;
mindt(1) = 35;
avdt(1) = 36.5;
maxdt(2) = 44.5;
mindt(2) = 41.5;
avdt(2) = 43;
maxdt(3) = 44.5;
mindt(3) = 41.5;
avdt(3) = 43;
maxdt(4) = 49;
mindt(4) = 46;
avdt(4) = 47.5;
maxdt(5) = 54.5;
mindt(5) = 50.5;
avdt(5) = 52.5;

ds = [1 2 3 5 10];

cmap = jet(11);

figure,
vwp = cell(5,1);
nlp = cell(5,1);
for i = 1 : 5
    subplot(1,5,i)
    hold on
    for j = 1 : length(str{i})
        errorbar(vw{i}(j), str{i}(j), errstr{i}(j),...
            'o', 'color', cmap(floor(nl{i}(j)/3)-2,:))
        if (str{i}(j)+errstr{i}(j)) <  maxdt(i) && ...
                (str{i}(j)-errstr{i}(j)) >  mindt(i)
            
            vwp{i} = vertcat(vwp{i}, vw{i}(j));
            nlp{i} = vertcat(nlp{i}, nl{i}(j));
        end
        
    end
    plot([0.1 25], [mindt(i) mindt(i)],'--g','linewidth',1)
    plot([0.1 25], [maxdt(i) maxdt(i)],'--g','linewidth',1)
    plot([0.1 25], [avdt(i) avdt(i)],'-g','linewidth',2)
    axis([0 25 10 75])
    title(['Dot Size: ' num2str(ds(i))])
end
%%
maxdt(1) = 38;
mindt(1) = 35;
avdt(1) = 36.5;
maxdt(2) = 44.5;
mindt(2) = 41.5;
avdt(2) = 43;
maxdt(3) = 44.5;
mindt(3) = 41.5;
avdt(3) = 43;
maxdt(4) = 49;
mindt(4) = 46;
avdt(4) = 47.5;
maxdt(5) = 54.5;
mindt(5) = 50.5;
avdt(5) = 52.5;

ds = [1 2 3 5 10];

cmap = jet(11);

figure,
vwp = cell(5,1);
nlp = cell(5,1);
for i = 1 : 5
    subplot(1,5,i)
    hold on
    for j = 1 : length(str{i})
%         if nl{i}(j) == 18
            errorbar(vw{i}(j), vi{i}(j), errvi{i}(j),...
                'o', 'color', cmap(floor(nl{i}(j)/3)-2,:))
            axis([0 25 -0.1 0.3])
%         end
    end
    title(['Dot Size: ' num2str(ds(i))])
end
%%
figure,
vwp = cell(5,1);
nlp = cell(5,1);
ds = [1 2 3 5 10];
for i = 1 : 5
%     subplot(1,5,i)
    hold on
    for j = 1 : length(str{i})
        if nl{i}(j) == 18 && vw{i}(j) == 6
            errorbar(ds(i), vi{i}(j), errvi{i}(j),...
                'o', 'color', cmap(floor(nl{i}(j)/3)-2,:))
            axis([0 11 -0.05 0.35])
        end
    end
    title(['Dot Size: ' num2str(ds(i))])
end



%%

figure, 
cmap = autumn(5);
hold on
for i = 1 : 5
    for j = 1 : length(str{i})
        if nl{i}(j) == 18 && vw{i}(j) == 6
            scatter(vi{i}(j), str{i}(j), 100, cmap(i,:), 'filled')
            
%             errorbar(ds(i), vi{i}(j), errvi{i}(j),...
%                 'o', 'color', cmap(nl{i}(j)/2-6,:))
            
            
            axis([-0.05 0.35 20 60])
        end
    end
    
end
%%
avdt(1) = 36.5;
avdt(2) = 41;
avdt(3) = 45;
avdt(4) = 47.5;
avdt(5) = 52.5;

figure, 
cmap = autumn(5);
hold on
estr = cell(5,1);
for i = 1 : 2
    for j = 1 : length(dt1.GMM(:,i))
            vi = (dt1.PSS2M(j,i)-dt1.PSS1M(j,i))/(dt1.PSS2M(j,i)+dt1.PSS1M(j,i));
            str = dt1.GMM(j,i);
            estr{i} = vertcat(estr{i}, str-avdt(i));
            scatter(vi, str, 100, cmap(i,:), 'filled')
            axis([-0.05 0.35 20 60])
    end
end
for i = 1 : 3
    for j = 1 : length(dt2.GMM(:,i))
            vi = (dt2.PSS2M(j,i)-dt2.PSS1M(j,i))/(dt2.PSS2M(j,i)+dt2.PSS1M(j,i));
            str = dt2.GMM(j,i);
            estr{i+2} = vertcat(estr{i+2}, str-avdt(i+2));
            scatter(vi, str, 100, cmap(2+i,:), 'filled')
            axis([-0.05 0.35 25 60])
    end
end

%%
figure,
hold on
cmap2 = jet(5);
for i = 1 : 5
    subplot(1,5,i)
    scatter(vwp{i}, nlp{i}, 100, cmap2(i,:), 'filled')
    axis([3 20 10 35])
end
%%
maxdt(1) = 0.03;
mindt(1) = 0.02;
avdt(1) = 0.025;
maxdt(2) = 0.095;
mindt(2) = 0.085;
avdt(2) = 0.09;
maxdt(3) = 0.095;
mindt(3) = 0.085;
avdt(3) = 0.09;
maxdt(4) = 0.23;
mindt(4) = 0.22;
avdt(4) = 0.225;
maxdt(5) = 0.295;
mindt(5) = 0.285;
avdt(5) = 0.29;

ds = [1 2 3 5 10];

cmap = jet(11);

figure,
for i = 1 : 5
    subplot(1,5,i)
    hold on
    for j = 1 : length(vi{i})
        errorbar(vw{i}(j), vi{i}(j), 0.05,...
            'o', 'color', cmap(nl{i}(j)/2-6,:))
    end
    plot([0.1 100], [mindt(i) mindt(i)],'--g','linewidth',1)
    plot([0.1 100], [maxdt(i) maxdt(i)],'--g','linewidth',1)
    plot([0.1 100], [avdt(i) avdt(i)],'-g','linewidth',2)
    axis([0 100 0 0.3])
    title(['Dot Size: ' num2str(ds(i))])
end

%%


%%
cmap = jet(11);
figure,
vwp = cell(5,1);
nlp = cell(5,1);
for i = 1 : 5
    subplot(1,5,i)
    hold on
    for j = 1 : length(str{i})
        errST = abs(str{i}(j)-avdt(i));
        
        
        
        errorbar(vw{i}(j), str{i}(j), errstr{i}(j),...
            'o', 'color', cmap(floor(nl{i}(j)/3)-2,:))
        if (str{i}(j)+errstr{i}(j)) <  maxdt(i) && ...
                (str{i}(j)-errstr{i}(j)) >  mindt(i)
            
            vwp{i} = vertcat(vwp{i}, vw{i}(j));
            nlp{i} = vertcat(nlp{i}, nl{i}(j));
        end
        
    end
    plot([0.1 100], [mindt(i) mindt(i)],'--g','linewidth',1)
    plot([0.1 100], [maxdt(i) maxdt(i)],'--g','linewidth',1)
    plot([0.1 100], [avdt(i) avdt(i)],'-g','linewidth',2)
    axis([0 25 10 75])
    title(['Dot Size: ' num2str(ds(i))])
end

%%
avdt(1) = 36.5;
avdt(2) = 41;
avdt(3) = 45;
avdt(4) = 47.5;
avdt(5) = 52.5;

ds = [1 2 3 5 10];
errVect = [];
for j = 1 : length(str{4})
    errST = 0;
    for i = 1 : 5
        errST = errST + (str{i}(j)-avdt(i))^2;
    end
    errST = sqrt(errST);
    errVect = vertcat(errVect,errST);
end
cmap = autumn(5);
cmap2 = winter(5);
ds = [1 2 3 5 10];
[B,I] = sort(errVect);

np = 16;
figure,
a = 1;
for j = 1 : length(nl{1})
%     if nl{1}(j) == 14
        subplot(4, ceil(np/4), a)
        a = a + 1;
        hold on
        for i = 1 : 5
            errorbar(ds(i), str{i}(j), errstr{i}(j), 'o', 'color', cmap(i,:))
            errorbar(ds(i), avdt(i), 1.5, 'o', 'color', cmap2(i,:))
            title(['nl ' num2str(nl{i}(j)) '; vw ' num2str(vw{i}(j)) ...
                '; mw ' num2str(mw{i}(j)) 'ERR ' num2str(errVect(j))])
        end
        axis([0 11 30 55])
%     end
end

% figure,
% np = 10;
% for j = 1 : np
%     subplot(2, ceil(np/2), j)
%     hold on
%     for i = 1 : 5
%         errorbar(ds(i), str{i}(I(j)), errstr{i}(I(j)), 'o', 'color', cmap(i,:))
%         errorbar(ds(i), avdt(i), 1.5, 'o', 'color', cmap2(i,:))
%         title(['nl ' num2str(nl{i}(I(j))) '; vw ' num2str(vw{i}(I(j))) ...
%             '; mw ' num2str(mw{i}(I(j))) 'ERR ' num2str(B(j))])
%     end
%     axis([0 11 30 55])
% end

%%
figure
hold on
in = find(B<10);
in = in(end);
for j = 1 : in
    for i = 1 : 5
        scatter(vw{i}(I(j)),nl{i}(I(j)),100,'ok')
    end
end

axis([0 100 12 35])
