path = 'C:\Users\Tomas\Desktop\Raw Low Res\VT058487 Cnt\';
flies = dir(path);
n = 3;
BoutsNG = cell(5,1);
BoutsRG = cell(5,1);
for n = 3 : length(flies)
    seq = load([path flies(n).name '\SaveProbs.mat'], 'seq');
    seq = seq.seq;
    dt = load([path flies(n).name '\DataLowRes.mat'],'Flies');
    dt = dt.Flies.Data;
    for k = 1 : length(seq)
        switch seq(k)
            case 1
                for i = 1 : length(dt{2*k-1}.Bouts);
                    vrb =  dt{2*k-1}.Vr(dt{2*k-1}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsNG{1});
                        BoutsNG{1}{l} = vrb;
                    end
                end
                for i = 1 : length(dt{2*k}.Bouts);
                    vrb =  dt{2*k}.Vr(dt{2*k}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsRG{1});
                        BoutsRG{1}{l} = vrb;
                    end
                end
            case 2.5
                for i = 1 : length(dt{2*k-1}.Bouts);
                    vrb =  dt{2*k-1}.Vr(dt{2*k-1}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsNG{2});
                        BoutsNG{2}{l} = vrb;
                    end
                end
                for i = 1 : length(dt{2*k}.Bouts);
                    vrb =  dt{2*k}.Vr(dt{2*k}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsRG{2});
                        BoutsRG{2}{l} = vrb;
                    end
                end
            case 5
                for i = 1 : length(dt{2*k-1}.Bouts);
                    vrb =  dt{2*k-1}.Vr(dt{2*k-1}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsNG{3});
                        BoutsNG{3}{l} = vrb;
                    end
                end
                for i = 1 : length(dt{2*k}.Bouts);
                    vrb =  dt{2*k}.Vr(dt{2*k}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsRG{3});
                        BoutsRG{3}{l} = vrb;
                    end
                end
            case 10
                for i = 1 : length(dt{2*k-1}.Bouts);
                    vrb =  dt{2*k-1}.Vr(dt{2*k-1}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsNG{4});
                        BoutsNG{4}{l} = vrb;
                    end
                end
                for i = 1 : length(dt{2*k}.Bouts);
                    vrb =  dt{2*k}.Vr(dt{2*k}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsRG{4});
                        BoutsRG{4}{l} = vrb;
                    end
                end
            case 15
                for i = 1 : length(dt{2*k-1}.Bouts);
                    vrb =  dt{2*k-1}.Vr(dt{2*k-1}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsNG{5});
                        BoutsNG{5}{l} = vrb;
                    end
                end
                for i = 1 : length(dt{2*k}.Bouts);
                    vrb =  dt{2*k}.Vr(dt{2*k}.Bouts{i});
                    if length(vrb) > 90
                        l = 1 + length(BoutsRG{5});
                        BoutsRG{5}{l} = vrb;
                    end
                end
                
            otherwise
                
        end
    end
    
    dd = cell(5,1);
    for i = 1 : 5
        for j = 1 : length(BoutsNG{i})
            %         visdns = BoutsNGVis{i}{j};
            visdns = BoutsNG{i}{j};
            vrdns = BoutsNG{i}{j}(1:length(visdns));
            dd{i} = vertcat(dd{i}, [vrdns visdns]);
            %         vrdns = downsample(BoutsNG{i}{j}, dsmp);
            %         visdns = downsample(BoutsNGVis{i}{j}, dsmp);
            %         dd{i} = vertcat(dd{i}, [vrdns visdns]);
        end
    end

    
end

%%
ang = [1 2 5 10 15];
nCopies = 10;
BoutsNGVis = cell(length(BoutsNG),1);
BoutsRGVis = cell(length(BoutsRG),1);
for j = 1 : length(BoutsNG)
    BoutsNGVis{j} = cell(length(BoutsNG{j}),1);
    for k = 1 : length(BoutsNG{j})
        vr = BoutsNG{j}{k};
        [Vresp, Vrp, Vmean, Vstd] = HRCPrediction(vr, ang(j), nCopies);
        BoutsNGVis{j}{k} = Vmean;
        disp(num2str(k));
    end
    for k = 1 : length(BoutsRG{j})
        vr = BoutsRG{j}{k};
        [Vresp, Vrp, Vmean, Vstd] = HRCPrediction(vr, ang(j), nCopies);
        BoutsRGVis{j}{k} = Vmean;
        disp(num2str(k));
    end
end

%%
dsmp = 15;
dd = cell(5,1);
for i = 1 : 5
    for j = 1 : length(BoutsNG{i})
%         visdns = BoutsNGVis{i}{j};
        visdns = BoutsNG{i}{j};
        vrdns = BoutsNG{i}{j}(1:length(visdns));
        dd{i} = vertcat(dd{i}, [vrdns visdns]);
%         vrdns = downsample(BoutsNG{i}{j}, dsmp);
%         visdns = downsample(BoutsNGVis{i}{j}, dsmp);
%         dd{i} = vertcat(dd{i}, [vrdns visdns]); 
    end
end

%
% figure,
% scatter(dd{1}(:,1), dd{1}(:,2), 20, 'b')
% hold on
% scatter(dd{2}(:,1), dd{2}(:,2), 20, 'c')
% hold on
% scatter(dd{3}(:,1), dd{3}(:,2), 20, 'g')
% hold on
% scatter(dd{4}(:,1), dd{4}(:,2), 20, 'r')
% hold on
% scatter(dd{5}(:,1), dd{5}(:,2), 20, 'k')
% 
% figure,
% scatter(dd{1}(2:end,2), dd{1}(1:end-1,1), 20, 'b')
% hold on
% scatter(dd{2}(2:end,2), dd{2}(1:end-1,1), 20, 'c')
% hold on
% scatter(dd{3}(2:end,2), dd{3}(1:end-1,1), 20, 'g')
% hold on
% scatter(dd{4}(2:end,2), dd{4}(1:end-1,1), 20, 'r')
% hold on
% scatter(dd{5}(2:end,2), dd{5}(1:end-1,1), 20, 'k')
%
vCents = 1000*(-1.5:0.05:1.5);
probsrg = nan*ones(length(vCents),5);
figure,
for j = 1 : 5
dd2 = [];
dd2(:,1) = dd{j}(16:end,1);
% dd2(:,2) = dd{j}(1:end-15,2);
dd2(:,2) = dd{j}(1:end-15,1);
% vCents = 0.001*(-2.5:0.05:2.5);

ns = zeros(length(vCents),1);
for k = 1 : length(vCents)-1
   aux1 = dd2(:,2);
   inds = find(aux1 > vCents(k) & aux1 <= vCents(k+1));
   aux2 = dd2(:,1);
   vals = aux2(inds);
   vals(vals > 0) = 1;
   vals(vals <= 0) = -1;
   probsrg(k,j) = sum(vals)/length(vals);
   ns(k) = length(vals);
end
l = floor(length(vCents)/2);

hold on
end
% axis([-0.002 0.002 -1 1])
plot(vCents, probsrg)
axis([-1000 1000 -1 1])

%%
dsmp = 15;
dd = cell(5,1);
for i = 1 : 5
   for j = 1 : length(BoutsRG{i})
%         vrdns = downsample(BoutsRG{i}{j}, dsmp);
%         visdns = downsample(BoutsRGVis{i}{j}, dsmp);
%         dd{i} = vertcat(dd{i}, [vrdns(1:length(visdns)) visdns]); 
%         visdns = BoutsRGVis{i}{j};
        visdns = BoutsRG{i}{j};
        vrdns = BoutsRG{i}{j}(1:length(visdns));
        dd{i} = vertcat(dd{i}, [vrdns visdns]);
   end
end

%
% figure,
% scatter(dd{1}(:,1), dd{1}(:,2), 20, 'b')
% hold on
% scatter(dd{2}(:,1), dd{2}(:,2), 20, 'c')
% hold on
% scatter(dd{3}(:,1), dd{3}(:,2), 20, 'g')
% hold on
% scatter(dd{4}(:,1), dd{4}(:,2), 20, 'r')
% hold on
% scatter(dd{5}(:,1), dd{5}(:,2), 20, 'k')
% title('Lag0')

% figure,
% scatter(dd{1}(2:end,2), dd{1}(1:end-1,1), 20, 'b')
% hold on
% scatter(dd{2}(2:end,2), dd{2}(1:end-1,1), 20, 'c')
% hold on
% scatter(dd{3}(2:end,2), dd{3}(1:end-1,1), 20, 'g')
% hold on
% scatter(dd{4}(2:end,2), dd{4}(1:end-1,1), 20, 'r')
% hold on
% scatter(dd{5}(2:end,2), dd{5}(1:end-1,1), 20, 'k')
% %%
%
vCents = 1000*(-1.5:0.05:1.5);
probsng = nan*ones(length(vCents),5);
figure,
for j = 1 : 5
dd2 = [];
dd2(:,1) = dd{j}(16:end,1);
% dd2(:,2) = dd{j}(1:end-15,2);
dd2(:,2) = dd{j}(1:end-15,1);
% vCents = 0.001*(-2.5:0.025:2.5);
ns = zeros(length(vCents),1);
for k = 1 : length(vCents)-1
   aux1 = dd2(:,2);
   inds = find(aux1 > vCents(k) & aux1 <= vCents(k+1));
   aux2 = dd2(:,1);
   vals = aux2(inds);
   vals(vals > 0) = 1;
   vals(vals <= 0) = -1;
   probsng(k,j) = sum(vals)/length(vals);
   ns(k) = length(vals);
end
l = floor(length(vCents)/2);
end
plot(vCents, probsng)
hold on
% axis([-0.002 0.002 -1 1])
axis([-1000 1000 -1 1])






