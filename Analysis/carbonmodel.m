clc;
clear all;
close all;
%mu = 0.01:0.005:0.05;
%beta = [0.03:0.01:0.09 0.1:0.05:0.25 0.3:0.1:0.5];
%alpha = 1.0:0.05:1.8;
%% now with some disturbances
fGain = @(AGB,a,b,c) a./(b+exp(-(AGB./c)));
fLoss = @(AGB,k) k.*AGB;

X     = 1000;
Y     = 1000;
N = 200;
arrsize = [X Y];
a   = 100.*ones(arrsize); % 100,b,500,0.1
%b  = max(normrnd(100,25,arrsize),1);
%b  = 0.1.*ones(arrsize);
e   = 500.*ones(arrsize);
k   = max(normrnd(0.05,0.05./4,arrsize),1E-3);
c0  = 10.*ones(arrsize);
mu    = 0.03;
alpha = 1.0;
I_slope = 0.5;

%% old cubes with wrong mu
% pwd = '/Net/Groups/BGI/people/swang/Data/D_tubes/Cubes_Intensity';  
% filename = sprintf('Cube_%dx%d_mu_%.3f_a_%.3f_Is_%.3f.mat',X,Y,mu,alpha,I_slope);
% cubedir = fullfile(pwd,filename);
% load(cubedir)

%% new cubes with right mu
in_cubedir = '/Net/Groups/BGI/scratch/swang/Data/D_cube/intensity/';
f_flag = sprintf('D_cube_Intenstiy_mu_%.3f_alpha_%.3f_beta_%.3f_on_*.mat',mu,alpha,I_slope);
f = dir(fullfile(in_cubedir,f_flag));
filepath = f.name;
cubedir = fullfile(in_cubedir,filepath);
load(cubedir)

D(D>=1)=1;
%figure,imagesc(D(:,:,5))

randIndex = randperm(size(D,3));
D_new = D(:,:,randIndex);

%% Display how different G change the GPP and AGB with disturbances
bs = [0.03,0.04,0.05,0.07,0.1];


figure(1)
%figure(2)
figure(2)
for i = 1:length(bs)
    b = bs(i);
    b_matrix = max(normrnd(b,b./4,arrsize),0.001);
    
    [c,g,n,l] = cdyn(c0,fGain,fLoss,a,b_matrix,e,k,N); % to see different Gs without disturbance
    %[c,g,n,l] = cdyn_DisT(c0,fGain,fLoss,a,b_matrix,e,k,N,D_new); 
    txt = ['G = ', num2str(b)];
    figure(1),plot(1:N,squeeze(mean(g,[1,2])),'DisplayName',txt,'LineWidth',2),hold on
    %title('GPP');
    ylabel('GPP (g\cdotm^{-2}\cdotyr^{-1})');
    xlabel('Year');
    %ylim([0,4000])
    % figure(2),plot(1:N,squeeze(mean(n,[1,2])),'DisplayName',txt,'LineWidth',2),hold on
    % %title('NPP');
    % ylabel('NPP (g\cdotm^{-2}\cdotyr^{-1})');
    % xlabel('Year');
    figure(2),plot(1:N,squeeze(mean(c,[1,2])),'DisplayName',txt,'LineWidth',2),hold on
    %title('AGB');
    ylabel('AGB (g\cdotm^{-2})');
    xlabel('Year');
    %ylim([0,40000])
    
end

%legend show
%close(2)


%% Display how different G change the GPP and AGB without disturbances
% figure(3)
% bs = [0.03,0.04,0.05,0.07,0.1];


% figure(4)
% figure(5)
% figure(6)
% for i = 1:length(bs)
%     b = bs(i);
%     [c,g,n,l] = cdyn(c0,fGain,fLoss,100,b,500,0.1,N); % (a b c)of Gain, k, N
%     txt = ['G = ', num2str(b)];
%     figure(4),plot(1:N,squeeze(mean(g,[1,2])),'DisplayName',txt),hold on
%     title('GPP')
%     figure(5),plot(1:N,squeeze(mean(n,[1,2])),'DisplayName',txt),hold on
%     title('NPP')
%     figure(6),plot(1:N,squeeze(mean(c,[1,2])),'DisplayName',txt),hold on
%     title('AGB')
    
% end

% legend show




% %%

% [cVeg_D,g,n,l] = cdyn_DisT(c0,fGain,fLoss,a,b,e,k,N,D);
% cVegFin = cVeg_D(:,:,end);
% %cVegFin = reshape(cVegFin,arrsize);



% %look the evolution of the cVeg map
% figure,colormap(jet)
% year     = [1,10,25,50,100,150];
% cVeg_max = prctile(cVeg_D(:),98);
% for ii = 1:size(year,2)
%     subplot(1,size(year,2),ii) %,hold on
%     cVeg_ = cVeg_D(:,:,year(ii));
%     imagesc(cVeg_)
%     title(['cVeg of year: ' num2str(year(ii))])
%     set(gca,'CLim',[0,cVeg_max]) % or caxis([cmin,cmax]);
%     axis image
% end
% colorbar('southoutside','Position',[0.1,0.3,0.8,0.05]); %[left, bottom, width, height],,'fontsize',15,

% curve = mean(cVeg_D,[1 2]);
% figure,plot(1:N,curve(:))

% cVeg_D = reshape(cVeg_D,[],N);
% figure, plot(1:N,cVeg_D(randperm(size(cVeg_D,1),3),:))
% ylabel('cVeg_D')
% %figure,plot(1:200,cVeg_D(randperm(size(cVeg_D,1),3),:))

% %the distribution of the final cVeg
% figure,histogram(cVegFin)


function [c,g,n,l] = cdyn(c0,fGain,fLoss,a,b,e,k,N)
    c       = NaN([size(c0),N]);
    g       = c;
    n       = c;
    l       = c;
    c(:,:,1)  = c0;
    for i = 2:N
        g(:,:,i-1)    = fGain(c(:,:,i-1),a,b,e);
        n(:,:,i-1)    = g(:,:,i-1).*0.5;
        l(:,:,i-1)	= fLoss(c(:,:,i-1),k);
        c(:,:,i)      = c(:,:,i-1)+n(:,:,i-1)-l(:,:,i-1);
    end
end

function [c,g,n,l] = cdyn_DisT(c0,fGain,fLoss,a,b,e,k,N,dcube) % a b e k should be matrix with size of size(c0)
    c           = NaN([size(c0),N]);
    g           = c;
    n           = c;
    l           = c;
    c(:,:,1)      = c0;
    
    
    for i = 2:N
        Cveg_init     = c(:,:,i-1).*(1-dcube(:,:,i-1));
        g(:,:,i-1)    = fGain(Cveg_init,a,b,e);
        n(:,:,i-1)    = g(:,:,i-1).*0.5;
        l(:,:,i-1)    = fLoss(Cveg_init,k);
        c(:,:,i)      = Cveg_init + n(:,:,i-1) - l(:,:,i-1); 
    end  

end
