function [pertVrMat, pertMat] = ProcessPerturbations(dt, PER, params)
PER = vertcat(abs(PER),0);
dpI = find(diff(PER)>0);
dpE = find(diff(PER)<0);

pertMat = [];
pertVrMat = [];
for i = 1 : length(dpI)
    if ((dpE(i)-dpI(i))>(params.pertD-3)) && (dpI(i)> 2*params.pertD) && ((dpI(i)+3*params.pertD) < length(dt.Vr))
        pertVrMat = horzcat(pertVrMat, dt.Vrus((dpI(i)-2*params.pertD):(dpI(i)+3*params.pertD)));
        pertMat = horzcat(pertMat, PER((dpI(i)-2*params.pertD):(dpI(i)+3*params.pertD)));
    end
end


% figure,
% hold on
% plot([-30 60], [0 0], '--g')
% plot([0 0], [-800 500], '--r')
% plot(-30:60, mean(pertVrMat,2)', 'k', 'linewidth', 2)
% plot(-30:60, mean(pertVrMat,2)' + std(pertVrMat,1,2)', 'k')
% plot(-30:60, mean(pertVrMat,2)' - std(pertVrMat,1,2)', 'k')
% plot(-30:60,-700+200*mean(pertMat,2)', 'b')
% plot(-30:60,-700+200*(mean(pertMat,2)' + std(pertMat,1,2)'), 'b')
% plot(-30:60,-700+200*(mean(pertMat,2)' - std(pertMat,1,2)'), 'b')
% axis([-30 60 -800 800])
% xlabel('Vr(º/s)')
% ylabel('Time (s)')
% title(['Perturbation :' params.PerP ' -> ' num2str(params.SizePerturbation)])


end