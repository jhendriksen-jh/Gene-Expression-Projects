%% Lab 2

%% using actual tumor data
clc, close all

r = zeros(50,1);
rlin = linspace(1,0,25);
r(1:25) = rlin;
g = zeros(50,1);
glin = linspace(0,1,26);
g(25:50) = glin;
map = zeros(50,3);
map(:,1) = r;
map(:,2) = g;

% I compiled a set of gene expression data from 35 female lung cancer
% samples using the TCGA with the data left in its original order

%load in data
tumordata = csvread('TCGA MicroArray Gene Expression Lung Female.csv',2,1);
tumordataSmoking = csvread('TCGA MicroArray Gene Expression Lung Female sorted by Smoking Indicator.csv',4,1);

figure(1)
%visualize data
subplot(1,2,1)
imagesc(tumordata)
subplot(1,2,2)
imagesc(tumordata)
colormap(map)


%% compute svd
[U,S,V] = svd(tumordata,0);
[Usmoke,Ssmoke,Vsmoke] = svd(tumordataSmoking,0);
V = V';
Vsmoke = Vsmoke';

imagesc(Ssmoke)
%%
%visualize SVD
figure(2)
subplot(1,2,1)
imagesc(U)
subplot(1,2,2)
imagesc(Usmoke)
colormap(map)

figure(3)
subplot(1,2,1)
imagesc(S)
subplot(1,2,2)
imagesc(Ssmoke)
colormap(map)

figure(4)
subplot(1,2,1)
imagesc(V)
subplot(1,2,2)
imagesc(Vsmoke)
colormap(map)

%% Visualize rows of V
figure(5)
subplot(1,2,1)
plot(V(4,:))
subplot(1,2,2)
plot(Vsmoke(4,:))

%% box plots for smoke of V' row 4
clc, close all
smokeind1 = Vsmoke(4,1:2);
smokeind2 = Vsmoke(4,3:6);
smokeind3 = Vsmoke(4,7:10);
smokeind4 = Vsmoke(4,11:35);
smoke = [smokeind1';smokeind2';smokeind3';smokeind4'];

i1 = repmat({'Non-Smoking'},length(smokeind1),1);
i2 = repmat({'Current Smoker'},length(smokeind2),1);
i3 = repmat({'Reformed Smoker > 15 years'},length(smokeind3),1);
i4 = repmat({'Reformed Smoker < 15 years'},length(smokeind4),1);
i = [i1;i2;i3;i4];

% boxplot(smoke,i)

pkw = kruskalwallis(smoke,i)
%% boxplots for current vs reformed smokers
smokeC = [smokeind2'];
smokeR = [smokeind3';smokeind4'];
Tc = repmat({'Current Smoker'},length(smokeind2),1);
Tr = repmat({'Reformed Smoker'},(length(smokeind3)+length(smokeind4)),1);
T = [Tc;Tr];

boxplot(smokeCvsR,T)
pCR = ranksum(smokeC,smokeR)
%% 
figure(6)
subplot(1,2,1)
imagesc(V(2:10,:))
subplot(1,2,2)
imagesc(Vsmoke(2:10,:))
colormap(map)


%% additional stuff 

% %% further processing of the data/svd
% %calc rank of S
% Srank = diag(Sreal);
% k=1;
% while k<=length(Srank)
%     if Srank(k,1) < Srank(1)/50
%         Srank(k,1) = 0;
%     end
%     k=k+1;
% end
% sigRank = nnz(Srank) 
% 
% %shape the data
% Ureal_ranked = Ureal(:,1:sigRank);
% Sreal_ranked = Sreal(1:sigRank,1:sigRank);
% Vreal_ranked = Vreal(:,1:sigRank);
% 
% %visualize shaped outputs
% figure(2)
% subplot(1,3,1)
% imagesc(Ureal)
% title('U')
% subplot(1,3,2)
% imagesc(Sreal)
% title('sigma')
% subplot(1,3,3)
% imagesc(Vreal')
% title('V')
% colormap(map)
% %%
% figure(3)
% subplot(2,3,1)
% imagesc(Ureal(:,2:9))
% title('U')
% subplot(2,3,2)
% imagesc(Sreal(2:9,2:9))
% title('sigma')
% subplot(2,3,3)
% imagesc(Vreal(:,2:9)')
% title('V')
% subplot(2,3,4)
% plot(Ureal(:,2:9))
% title('U')
% subplot(2,3,5)
% plot(Sreal(2:9,2:9))
% title('sigma')
% subplot(2,3,6)
% plot(Vreal(:,2:9))
% title('V')
% colormap(map)
% %% sort v?
% vsort = sort(Vreal);
% rowVsort = sortrows(Vreal,2);
% subplot(1,3,1)
% imagesc(Vreal')
% subplot(1,3,2)
% imagesc(vsort')
% subplot(1,3,3)
% imagesc(rowVsort')
% colormap(map)
% 
% figure(4)
% subplot(3,1,1)
% plot(vsort(:,1:5))
% subplot(3,1,2)
% plot(vsort(:,6:10))
% subplot(3,1,3)
% plot(vsort(:,11:15))
% %% testing svd outputs
% %multiply outputs to generate input
% testfull = Ureal*Sreal*Vreal';
% testrank = Ureal_ranked*Sreal_ranked*Vreal_ranked';
% %% testing orthonormality
% clc, close all
% 
% Utest = Ureal'*Ureal;
% Vtest = Vreal*Vreal';
% subplot(1,2,1)
% imshow(Utest)
% title('U')
% subplot(1,2,2)
% imshow(Vtest)
% title('V')
% 
% %% Visualizing components of the SVD from the real data
% clc, close all
% %visualizing tests of svd outputs
% figure(3)
% subplot(3,1,1)
% imagesc(tumordata)
% title('input data')
% subplot(3,1,2)
% imagesc(testfull)
% title('test of full outputs')
% subplot(3,1,3)
% imagesc(testrank)
% title('test of shaped outputs')
% colormap(map)
% 
% %visualzing the first 3 signals in V and U
% figure(5)
% subplot(3,1,1)
% plot(Vreal_ranked(:,1))
% title('V1')
% subplot(3,1,2)
% plot(Vreal_ranked(:,2))
% title('V2')
% subplot(3,1,3)
% plot(Vreal_ranked(:,3))
% title('V3')
% 
% figure(6)
% subplot(3,1,1)
% plot(Ureal_ranked(:,1))
% title('U1')
% subplot(3,1,2)
% plot(Ureal_ranked(:,2))
% title('U2')
% subplot(3,1,3)
% plot(Ureal_ranked(:,3))
% title('U3')