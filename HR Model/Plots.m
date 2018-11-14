ang = 1;
Vr = Dt.VrNG{1}{1};
NCopies = 10;
[Vresp1, Vrp] = HRCPrediction(Vr, ang, NCopies);
ang = 2;
% Vr = Dt.VrNG{2}{1};
% NCopies = 10;
[Vresp2, Vrp] = HRCPrediction(Vr, ang, NCopies);
ang = 5;
% Vr = Dt.VrNG{3}{1};
% NCopies = 10;
[Vresp5, Vrp] = HRCPrediction(Vr, ang, NCopies);
ang = 10;
% Vr = Dt.VrNG{4}{1};
% NCopies = 10;
[Vresp10, Vrp] = HRCPrediction(Vr, ang, NCopies);
% ang = 15;
% Vr = Dt.VrNG{5}{1};
% NCopies = 10;
% [Vresp15, Vrp] = HRCPrediction(Vr, ang, NCopies);

%
c = jet(10);
figure,
plot(mean(Vresp1,2),'color', c(1,:), 'linewidth',2)
hold on
plot(mean(Vresp2,2),'color', c(3,:), 'linewidth',2)
hold on
plot(mean(Vresp5,2),'color', c(7,:), 'linewidth',2)
hold on
plot(mean(Vresp10,2),'color', c(9,:), 'linewidth',2)
hold on
% plot(mean(Vresp15,2),'color', c(9,:), 'linewidth',2)
% hold on

%%
vrCents = -20:1:20;
m = mean(Vresp1,2);
m = m(1:end-1);
vMean1 = [];
vStd1 = [];
for j = 1 : length(vrCents)-1;
    indx = find(Vrp==vrCents(j));
    if length(indx) > 5
        vMean1 = vertcat(vMean1, mean(m(indx)));
        vStd1 = vertcat(vStd1, std(m(indx)));
    else
        vMean1 = vertcat(vMean1, nan);
        vStd1 = vertcat(vStd1, nan);
    end
end
m = mean(Vresp2,2);
m = m(1:end-1);
vMean2 = [];
vStd2 = [];
for j = 1 : length(vrCents)-1;
    indx = find(Vrp==vrCents(j));
    if length(indx) > 5
        vMean2 = vertcat(vMean2, mean(m(indx)));
        vStd2 = vertcat(vStd2, std(m(indx)));
    else
        vMean2 = vertcat(vMean2, nan);
        vStd2 = vertcat(vStd2, nan);
    end
end
m = mean(Vresp5,2);
m = m(1:end-1);
vMean5 = [];
vStd5 = [];
for j = 1 : length(vrCents)-1;
    indx = find(Vrp==vrCents(j));
    if length(indx) > 5
        vMean5 = vertcat(vMean5, mean(m(indx)));
        vStd5 = vertcat(vStd5, std(m(indx)));
    else
        vMean5 = vertcat(vMean5, nan);
        vStd5 = vertcat(vStd5, nan);
    end
end
m = mean(Vresp10,2);
m = m(1:end-1);
vMean10 = [];
vStd10 = [];
for j = 1 : length(vrCents)-1;
    indx = find(Vrp==vrCents(j));
    if length(indx) > 5
        vMean10 = vertcat(vMean10, mean(m(indx)));
        vStd10 = vertcat(vStd10, std(m(indx)));
    else
        vMean10 = vertcat(vMean10, nan);
        vStd10 = vertcat(vStd10, nan);
    end
end
% m = mean(Vresp15,2);
% m = m(1:end-1);
% vMean15 = [];
% vStd15 = [];
% for j = 1 : length(vrCents)-1;
%     indx = find(Vrp==vrCents(j));
%     if length(indx) > 5
%         vMean15 = vertcat(vMean15, mean(m(indx)));
%         vStd15 = vertcat(vStd15, std(m(indx)));
%     else
%         vMean15 = vertcat(vMean15, nan);
%         vStd15 = vertcat(vStd15, nan);
%     end
% end

figure,
errorbar(vrCents(2:end), vMean1, vStd1, 'color', c(1,:), 'linewidth', 2)
hold on
errorbar(vrCents(2:end), vMean2, vStd2, 'color', c(3,:), 'linewidth', 2)
hold on
errorbar(vrCents(2:end), vMean5, vStd5, 'color', c(5,:), 'linewidth', 2)
hold on
errorbar(vrCents(2:end), vMean10, vStd10, 'color', c(7,:), 'linewidth', 2)
hold on
% errorbar(vrCents(2:end), vMean15, vStd15, 'color', c(10,:), 'linewidth', 2)
% hold on
%%
dsmp = 15;
Vrdns = downsample(Vr, dsmp);
m = mean(Vresp5,2);
m = m(1:end-1);
mdns = downsample(m, dsmp);
figure,
plot([min(mdns(1:end-1)) max(mdns(1:end-1))], zeros(2,1), 'g', 'linewidth', 2)
hold on
plot(zeros(2,1), [min(Vrdns(1:end-1)) max(Vrdns(1:end-1))], 'r', 'linewidth', 2)
hold on
scatter(mdns(1:end-1),Vrdns(2:end),100, 'k')