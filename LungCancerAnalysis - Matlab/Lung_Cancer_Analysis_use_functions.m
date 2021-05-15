%%
clc, close all
[pvmww,pkv] = svdBoxPlotsSmoke('TCGA MicroArray Gene Expression Lung Female sorted by Smoking Indicator.csv',2,35,0.2);
pkv
%%
clc, close all, clear all

Cat = [{'Non-Smoking'},{'Current Smoker'},{'Reformed > 15 yrs'},{'Reformed < 15 yrs'}];
CatL = [3,4,4,25];

[pvbp,pkv] = svdBoxPlotsGen('TCGA MicroArray Gene Expression Lung Female sorted by Smoking Indicator.csv',4,35,0.05,Cat,CatL);
pkv
%% Visualize Women vs Combined Datasets (combined have less rows)
r = zeros(50,1);
rlin = linspace(0,1,25);
r(26:50) = rlin;
g = zeros(50,1);
glin = linspace(1,0,25);
g(1:25) = glin;
map = zeros(50,3);
map(:,1) = r;
map(:,2) = g;

%load in data
Wdata = csvread('TCGA MicroArray Gene Expression Lung Female.csv',2,1);
Cdata = csvread('Combined Lung Genes.csv',2,1);

figure(1)
%visualize data
subplot(1,2,1)
imagesc(Wdata)
subplot(1,2,2)
imagesc(Cdata)
colormap(map)
%%
imagesc(Cdata)
colormap(map)
title('Rasterization of Gene Expression Data')
xlabel('Patient')
ylabel('Gene')

%%
[Uw,Sw,Vw] = svd(Wdata,0);
[Uc,Sc,Vc] = svd(Cdata,0);
Vw = Vw';
Vc = Vc';
%%
figure(1)
Uc2 = Uc(:,1:130);
imagesc(Uc2)
title("U")
colormap(map)

figure(2)
plot(diag(Sc),'-o')
title("Sigma")
colormap(map)

figure(3)
imagesc(Vc)
title("V transpose")
colormap(map)
%%
figure(2)
subplot(2,3,1)
imagesc(Uw)
subplot(2,3,2)
imagesc(Sw)
subplot(2,3,3)
imagesc(Vw)
subplot(2,3,4)
imagesc(Uc)
subplot(2,3,5)
imagesc(Sc)
subplot(2,3,6)
imagesc(Vc)
colormap(map)

figure(3)
imagesc(Vc)
colormap(map)
figure(4)
imagesc(Vw)
colormap(map)
%% barchart sigma
figure(5)
bar(Sw,20)
figure(6)
bar(Sc,80)
Scdiag = diag(Sc)
Swdiag = diag(Sw)
%% Combined vs Combined Sorted
r = zeros(50,1);
rlin = linspace(0,1,25);
r(26:50) = rlin;
g = zeros(50,1);
glin = linspace(1,0,25);
g(1:25) = glin;
map = zeros(50,3);
map(:,1) = r;
map(:,2) = g;

%load in data
Sdata = csvread('CombinedSortedSmoking.csv',4,1);
Cdata = csvread('Combined Lung Genes.csv',2,1);
Sexdata = csvread('CombinedSortedSex.csv',4,1);
PSdata = csvread('CombinedSortedPrimaryStage.csv',4,1);
Mdata = csvread('CombinedSortedMetastasis.csv',4,1);
Hdata = csvread('CombinedSortedHistory.csv',4,1);
Adata = csvread('CombinedSortedAge.csv',4,1);

%%
figure(1)
%visualize data
subplot(1,2,1)
imagesc(Sdata)
subplot(1,2,2)
imagesc(Cdata)
colormap(map)

[Us,Ss,Vs] = svd(Sdata,0);
[Uc,Sc,Vc] = svd(Cdata,0);
[Usex,Ssex,Vsex] = svd(Sexdata,0);
[Ups,Sps,Vps] = svd(PSdata,0);
[Um,Sm,Vm] = svd(Mdata,0);
[Uh,Sh,Vh] = svd(Hdata,0);
[Ua,Sa,Va] = svd(Adata,0);

Vs = Vs';
Vc = Vc';
Vsex = Vsex';
Vps = Vps';
Vm = Vm';
Vh = Vh';
Va = Va';
%%
figure(1)
Usex2 = Usex(:,1:130);
imagesc(Usex2)
title("U")
colormap(map)

figure(2)
plot(diag(Ssex),'-o')
title("Sigma")
colormap(map)

figure(3)
imagesc(Vsex)
title("V transpose")
colormap(map)
%%
figure(1)
imagesc(Vc(1:5,:))
colormap(map)
title('unsorted')
figure(2)
imagesc(Vsex(1:5,:))
colormap(map)
title('sorted by gender')
%% exporting U's as csv's
csvwrite('Usex.csv',Usex)
%%
figure(2)

subplot(2,3,1)
imagesc(Us)
subplot(2,3,2)
imagesc(Ss)
subplot(2,3,3)
imagesc(Vs)
subplot(2,3,4)
imagesc(Uc)
subplot(2,3,5)
imagesc(Sc)
subplot(2,3,6)
imagesc(Vc)
colormap(map)
%%
figure(3)
subplot(2,3,1)
imagesc(Vs(1:20,:))
title('smoking')
subplot(2,3,2)
imagesc(Vsex(1:20,:))
title('sex')
subplot(2,3,3)
imagesc(Vps(1:20,:))
title('Stage Tumor')
subplot(2,3,4)
imagesc(Vm(1:20,:))
title('Metastasis')
subplot(2,3,5)
imagesc(Vh(1:20,:))
title('History')
subplot(2,3,6)
imagesc(Va(1:20,:))
title('age')
colormap(map)
%%
clc, close all, clear all

Cat = [{'Non-Smoking'},{'Current Smoker'},{'Reformed > 15 yrs'},{'Reformed < 15 yrs'}];
CatL = [5,20,27,80];

[pvbp,pkv] = svdBoxPlotsSmoke('CombinedSortedSmoking.csv',4,132,0.05,Cat,CatL);
pkv

%%
clc, close all,clear all
Cat = [{'Non-Smoking'},{'Current or Former Smoker'}];
CatL = [5,127];

[pvbp,pkv] = svdBoxPlotsGen('CombinedSortedSmoking.csv',2,132,0.1,Cat,CatL);
pkv
%% Using function on all combined sorted


% sorted by smoking
CatS = [{'Non-Smoking'},{'Current or Former Smoker'}];
LS = [5,127];
[pVkS,pVS] = svdBoxPlotsGen('CombinedSortedSmoking.csv',2,132,0.01,CatS,LS);
%%
% sorted by Sex
CatSex = [{'Female'},{'Male'}];
LSex = [35,97];
[pVkSex,pVSex] = svdBoxPlotsGen('CombinedSortedSex.csv',2,132,0.01,CatSex,LSex);
pVSex(2)
DS = diag(Ssex);
S2 = DS(1:3)
%%
% sorted by Tumor Stage
CatPS = [{'T1'},{'T2'},{'T3'},{'T4'}];
LPS = [22,87,12,11];
[pVkPS,pVPS] = svdBoxPlotsGen('CombinedSortedPrimaryStage.csv',4,132,0.01,CatPS,LPS);

CatPS2 = [{'Low Stage'},{'High Stage'}];
LPS2 = [22,110];
[pVkPS2,pVPS2] = svdBoxPlotsGen('CombinedSortedPrimaryStage.csv',2,132,0.01,CatPS2,LPS2);

pVPS2
%%
DS = diag(Sps);
%%
csvwrite('UEarlyLateStage.csv',Ups);
%%
csvwrite('USmoking.csv',Us);
csvwrite('UageDiag.csv',Ua);
%%
% sorted by Metastasis
CatM = [{'No Metastasis'},{'Metastasis'}];
LM = [129,3];
[pVkM,pVM] = svdBoxPlotsGen('CombinedSortedMetastasis.csv',2,132,0.05,CatM,LM);

% sorted by History
CatH = [{'No History'},{'History'}];
LH = [109,23];
[pVkH,pVH] = svdBoxPlotsGen('CombinedSortedHistory.csv',2,132,0.05,CatH,LH);
%%
% sorted by age
CatA = [{'Over 67'},{'Under 67'}];
LA = [64,68];
[pVkA,pVA] = svdBoxPlotsGen('CombinedSortedAge.csv',2,132,0.05,CatA,LA);
pVA
%%
% sorted by age with 4 categories
CatAex = [{'Over 70'},{'66-70'},{'60-65'},{'Under 60'}];
LAex = [46,24,36,26];
[pVkAex,pVAex] = svdBoxPlotsGen('CombinedSortedAge.csv',3,132,0.05,CatAex,LAex);

%% 
close all, clc
pVS
S8 = DS(7:9)
S15 = DS(14:16)
pVSex
pVPS
PS4 = DS(3:5)
pVM
M31 = DS(30:32)
pVH
H22 = DS(21:23)
pVA
A7 = DS(6:8)
A12 = DS(11:13)
A23 = DS(22:24)

%%
t2500 = hygecdf(31,12024,123,2500)
t100 = hygecdf(13,12024,123,1000)
t2500 = hygecdf(3,12024,123,250)



