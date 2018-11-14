d1 = load('C:\Users\tomas\Desktop\SimAB8.51.mat');
d5 = load('C:\Users\tomas\Desktop\SimAB8.55.mat');
d10 = load('C:\Users\tomas\Desktop\SimAB8.510.mat');
gMNG = zeros(3,1);
gMRG = zeros(3,1);
semMNG = zeros(3,1);
semMRG = zeros(3,1);

gMNG(1) = (d1.dtt.gmvFTBNG-d1.dtt.gmvBTFNG)/(d1.dtt.gmvBTFNG+d1.dtt.gmvFTBNG);
gMNG(2) = (d5.dtt.gmvFTBNG-d5.dtt.gmvBTFNG)/(d5.dtt.gmvBTFNG+d5.dtt.gmvFTBNG);
gMNG(3) = (d10.dtt.gmvFTBNG-d10.dtt.gmvBTFNG)/(d10.dtt.gmvBTFNG+d10.dtt.gmvFTBNG);

semMNG(1) = sqrt(((2*d1.dtt.gmvBTFNG/((d1.dtt.gmvBTFNG+d1.dtt.gmvFTBNG)^2))*d1.dtt.semvFTBNG)^2 + ...
    ((2*d1.dtt.gmvFTBNG/((d1.dtt.gmvBTFNG+d1.dtt.gmvFTBNG)^2))*d1.dtt.semvBTFNG)^2);
semMNG(2) = sqrt(((2*d5.dtt.gmvBTFNG/((d5.dtt.gmvBTFNG+d5.dtt.gmvFTBNG)^2))*d5.dtt.semvFTBNG)^2 + ...
    ((2*d5.dtt.gmvFTBNG/((d5.dtt.gmvBTFNG+d5.dtt.gmvFTBNG)^2))*d5.dtt.semvBTFNG)^2);
semMNG(3) = sqrt(((2*d10.dtt.gmvBTFNG/((d10.dtt.gmvBTFNG+d10.dtt.gmvFTBNG)^2))*d10.dtt.semvFTBNG)^2 + ...
    ((2*d10.dtt.gmvFTBNG/((d10.dtt.gmvBTFNG+d10.dtt.gmvFTBNG)^2))*d10.dtt.semvBTFNG)^2);

gMRG(1) = (d1.dtt.gmvFTBRG-d1.dtt.gmvBTFRG)/(d1.dtt.gmvBTFRG+d1.dtt.gmvFTBRG);
gMRG(2) = (d5.dtt.gmvFTBRG-d5.dtt.gmvBTFRG)/(d5.dtt.gmvBTFRG+d5.dtt.gmvFTBRG);
gMRG(3) = (d10.dtt.gmvFTBRG-d10.dtt.gmvBTFRG)/(d10.dtt.gmvBTFRG+d10.dtt.gmvFTBRG);

semMRG(1) = sqrt(((2*d1.dtt.gmvBTFRG/((d1.dtt.gmvBTFRG+d1.dtt.gmvFTBRG)^2))*d1.dtt.semvFTBRG)^2 + ...
    ((2*d1.dtt.gmvFTBRG/((d1.dtt.gmvBTFRG+d1.dtt.gmvFTBRG)^2))*d1.dtt.semvBTFRG)^2);
semMRG(2) = sqrt(((2*d5.dtt.gmvBTFRG/((d5.dtt.gmvBTFNG+d5.dtt.gmvFTBRG)^2))*d5.dtt.semvFTBRG)^2 + ...
    ((2*d5.dtt.gmvFTBRG/((d5.dtt.gmvBTFRG+d5.dtt.gmvFTBRG)^2))*d5.dtt.semvBTFRG)^2);
semMRG(3) = sqrt(((2*d10.dtt.gmvBTFRG/((d10.dtt.gmvBTFRG+d10.dtt.gmvFTBRG)^2))*d10.dtt.semvFTBRG)^2 + ...
    ((2*d10.dtt.gmvFTBRG/((d10.dtt.gmvBTFRG+d10.dtt.gmvFTBRG)^2))*d10.dtt.semvBTFRG)^2);

figure,
errorbar([1 5 10], gMNG, semMNG, 'o', 'color', 'r', 'linewidth', 3)
hold on
errorbar([1 5 10], gMRG, semMRG, 'o', 'color', 'b', 'linewidth', 3)

axis([0 12 -0.1 0.3])
xlabel('Dot Size [º]')
ylabel('Direction Asymmetry [a.u.]')




%%

d1 = load('C:\Users\tomas\Desktop\SimAB12.36_21.mat');
d5 = load('C:\Users\tomas\Desktop\SimAB12.36_25.mat');
d10 = load('C:\Users\tomas\Desktop\SimAB12.36_210.mat');
gMNG = zeros(3,1);
gMRG = zeros(3,1);
semMNG = zeros(3,1);
semMRG = zeros(3,1);

gMNG(1) = (d1.dtt.gmvFTBNG-d1.dtt.gmvBTFNG)/(d1.dtt.gmvBTFNG+d1.dtt.gmvFTBNG);
gMNG(2) = (d5.dtt.gmvFTBNG-d5.dtt.gmvBTFNG)/(d5.dtt.gmvBTFNG+d5.dtt.gmvFTBNG);
gMNG(3) = (d10.dtt.gmvFTBNG-d10.dtt.gmvBTFNG)/(d10.dtt.gmvBTFNG+d10.dtt.gmvFTBNG);

semMNG(1) = sqrt(((2*d1.dtt.gmvBTFNG/((d1.dtt.gmvBTFNG+d1.dtt.gmvFTBNG)^2))*d1.dtt.semvFTBNG)^2 + ...
    ((2*d1.dtt.gmvFTBNG/((d1.dtt.gmvBTFNG+d1.dtt.gmvFTBNG)^2))*d1.dtt.semvBTFNG)^2);
semMNG(2) = sqrt(((2*d5.dtt.gmvBTFNG/((d5.dtt.gmvBTFNG+d5.dtt.gmvFTBNG)^2))*d5.dtt.semvFTBNG)^2 + ...
    ((2*d5.dtt.gmvFTBNG/((d5.dtt.gmvBTFNG+d5.dtt.gmvFTBNG)^2))*d5.dtt.semvBTFNG)^2);
semMNG(3) = sqrt(((2*d10.dtt.gmvBTFNG/((d10.dtt.gmvBTFNG+d10.dtt.gmvFTBNG)^2))*d10.dtt.semvFTBNG)^2 + ...
    ((2*d10.dtt.gmvFTBNG/((d10.dtt.gmvBTFNG+d10.dtt.gmvFTBNG)^2))*d10.dtt.semvBTFNG)^2);

gMRG(1) = (d1.dtt.gmvFTBRG-d1.dtt.gmvBTFRG)/(d1.dtt.gmvBTFRG+d1.dtt.gmvFTBRG);
gMRG(2) = (d5.dtt.gmvFTBRG-d5.dtt.gmvBTFRG)/(d5.dtt.gmvBTFRG+d5.dtt.gmvFTBRG);
gMRG(3) = (d10.dtt.gmvFTBRG-d10.dtt.gmvBTFRG)/(d10.dtt.gmvBTFRG+d10.dtt.gmvFTBRG);

semMRG(1) = sqrt(((2*d1.dtt.gmvBTFRG/((d1.dtt.gmvBTFRG+d1.dtt.gmvFTBRG)^2))*d1.dtt.semvFTBRG)^2 + ...
    ((2*d1.dtt.gmvFTBRG/((d1.dtt.gmvBTFRG+d1.dtt.gmvFTBRG)^2))*d1.dtt.semvBTFRG)^2);
semMRG(2) = sqrt(((2*d5.dtt.gmvBTFRG/((d5.dtt.gmvBTFNG+d5.dtt.gmvFTBRG)^2))*d5.dtt.semvFTBRG)^2 + ...
    ((2*d5.dtt.gmvFTBRG/((d5.dtt.gmvBTFRG+d5.dtt.gmvFTBRG)^2))*d5.dtt.semvBTFRG)^2);
semMRG(3) = sqrt(((2*d10.dtt.gmvBTFRG/((d10.dtt.gmvBTFRG+d10.dtt.gmvFTBRG)^2))*d10.dtt.semvFTBRG)^2 + ...
    ((2*d10.dtt.gmvFTBRG/((d10.dtt.gmvBTFRG+d10.dtt.gmvFTBRG)^2))*d10.dtt.semvBTFRG)^2);

figure,
errorbar([1 5 10], gMNG, semMNG, 'o', 'color', 'r', 'linewidth', 3)
hold on
errorbar([1 5 10], gMRG, semMRG, 'o', 'color', 'b', 'linewidth', 3)

axis([0 12 -0.1 0.3])
xlabel('Dot Size [º]')
ylabel('Direction Asymmetry [a.u.]')

%%

d1 = load('C:\Users\tomas\Desktop\SimAB8.51.mat');
d5 = load('C:\Users\tomas\Desktop\SimAB8.55.mat');
d10 = load('C:\Users\tomas\Desktop\SimAB8.510.mat');
gMNG = zeros(3,1);
gMRG = zeros(3,1);
semMNG = zeros(3,1);
semMRG = zeros(3,1);

gMNG(1) = (d1.dtt.gmvFTBNG-d1.dtt.gmvBTFNG)/(d1.dtt.gmvBTFNG+d1.dtt.gmvFTBNG);
gMNG(2) = (d5.dtt.gmvFTBNG-d5.dtt.gmvBTFNG)/(d5.dtt.gmvBTFNG+d5.dtt.gmvFTBNG);
gMNG(3) = (d10.dtt.gmvFTBNG-d10.dtt.gmvBTFNG)/(d10.dtt.gmvBTFNG+d10.dtt.gmvFTBNG);

semMNG(1) = sqrt(((2*d1.dtt.gmvBTFNG/((d1.dtt.gmvBTFNG+d1.dtt.gmvFTBNG)^2))*d1.dtt.semvFTBNG)^2 + ...
    ((2*d1.dtt.gmvFTBNG/((d1.dtt.gmvBTFNG+d1.dtt.gmvFTBNG)^2))*d1.dtt.semvBTFNG)^2);
semMNG(2) = sqrt(((2*d5.dtt.gmvBTFNG/((d5.dtt.gmvBTFNG+d5.dtt.gmvFTBNG)^2))*d5.dtt.semvFTBNG)^2 + ...
    ((2*d5.dtt.gmvFTBNG/((d5.dtt.gmvBTFNG+d5.dtt.gmvFTBNG)^2))*d5.dtt.semvBTFNG)^2);
semMNG(3) = sqrt(((2*d10.dtt.gmvBTFNG/((d10.dtt.gmvBTFNG+d10.dtt.gmvFTBNG)^2))*d10.dtt.semvFTBNG)^2 + ...
    ((2*d10.dtt.gmvFTBNG/((d10.dtt.gmvBTFNG+d10.dtt.gmvFTBNG)^2))*d10.dtt.semvBTFNG)^2);

gMRG(1) = (d1.dtt.gmvFTBRG-d1.dtt.gmvBTFRG)/(d1.dtt.gmvBTFRG+d1.dtt.gmvFTBRG);
gMRG(2) = (d5.dtt.gmvFTBRG-d5.dtt.gmvBTFRG)/(d5.dtt.gmvBTFRG+d5.dtt.gmvFTBRG);
gMRG(3) = (d10.dtt.gmvFTBRG-d10.dtt.gmvBTFRG)/(d10.dtt.gmvBTFRG+d10.dtt.gmvFTBRG);

semMRG(1) = sqrt(((2*d1.dtt.gmvBTFRG/((d1.dtt.gmvBTFRG+d1.dtt.gmvFTBRG)^2))*d1.dtt.semvFTBRG)^2 + ...
    ((2*d1.dtt.gmvFTBRG/((d1.dtt.gmvBTFRG+d1.dtt.gmvFTBRG)^2))*d1.dtt.semvBTFRG)^2);
semMRG(2) = sqrt(((2*d5.dtt.gmvBTFRG/((d5.dtt.gmvBTFNG+d5.dtt.gmvFTBRG)^2))*d5.dtt.semvFTBRG)^2 + ...
    ((2*d5.dtt.gmvFTBRG/((d5.dtt.gmvBTFRG+d5.dtt.gmvFTBRG)^2))*d5.dtt.semvBTFRG)^2);
semMRG(3) = sqrt(((2*d10.dtt.gmvBTFRG/((d10.dtt.gmvBTFRG+d10.dtt.gmvFTBRG)^2))*d10.dtt.semvFTBRG)^2 + ...
    ((2*d10.dtt.gmvFTBRG/((d10.dtt.gmvBTFRG+d10.dtt.gmvFTBRG)^2))*d10.dtt.semvBTFRG)^2);

figure,
errorbar([1 5 10], gMNG, semMNG, 'o', 'color', 'r', 'linewidth', 3)
hold on
errorbar([1 5 10], gMRG, semMRG, 'o', 'color', 'b', 'linewidth', 3)

axis([0 12 -0.1 0.3])
xlabel('Dot Size [º]')
ylabel('Direction Asymmetry [a.u.]')
