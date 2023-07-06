
% DisturbanceRegimeParameters : mu, alpha, Iz(currently only slope);
% Biomass Statistics: mean, median, range, skewness, kurtosis.
function Calculate_Stats(mu,alpha,beta,GPP_index,Lb)
    %mu= [0.01:0.005:0.05];
    %alpha = [1.0:0.05:1.8];
    %beta = [0.03:0.01:0.09 0.1:0.05:0.25 0.3:0.1:0.5];
    %GPP_index = [0.13 0.14 0.15]
     
    time1 = datetime('now');

    fGain = @(AGB,a,b,c) a./(b+exp(-(AGB./c)));
    fLoss = @(AGB,k) k.*AGB;

    X     = 1000;
    Y     = 1000;
    N     = 200;
    arrsize = [X Y];
    a   = 100.*ones(arrsize);
    %a   = max(normrnd(100,25,arrsize),1); % 100,b,500,0.1
    b  = max(normrnd(GPP_index,GPP_index./4,arrsize),1E-3);  %GPP index:[0.05:0.01:0.15]
    %b   = 0.1.*ones(arrsize);
    e   = 500.*ones(arrsize);%fixed scale
    k   = max(normrnd(Lb, Lb./4,arrsize),1E-3);
    c0  = 10.*ones(arrsize);

    %outdirs and outnames 
    %cvegdir = '/Net/Groups/BGI/scratch/swang/Data/Paper2/Statistics/cveg/';
    %gppdir = '/Net/Groups/BGI/scratch/swang/Data/Paper2/Statistics/gpp/'; 
    xlsdir = sprintf('./xls/k_%.3f/', Lb);
    
    if ~exist(xlsdir,'dir')
        mkdir(xlsdir)
    end

    %cveg
    % out_cvegname = sprintf('cVegFin_%dx%d_Mu_%.3f_Alpha_%.3f_Beta_%.3f_G_%.2f.mat',X,Y,mu,alpha,beta,GPP_index);
    % out_cvegfullname = fullfile(cvegdir,out_cvegname);
    % %gpp
    % out_gppname = sprintf('GPP_%dx%d_Mu_%.3f_Alpha_%.3f_Beta_%.3f_G_%.2f.mat',X,Y,mu,alpha,beta,GPP_index);
    % out_gppfullname = fullfile(gppdir,out_gppname);
    %xls
    xlsname = sprintf('xls_%dx%d_Mu_%.3f_Alpha_%.3f_Beta_%.3f_G_%.2f_K_%.2f_all_statistics.xls',X,Y,mu,alpha,beta,GPP_index,Lb);
    xlsfile = fullfile(xlsdir,xlsname);

    if exist(xlsfile, 'file') == 2
        disp(xlsfile)
    else
        T_mu = [];
        T_Is = []; 
        T_alpha = [];
        T_kurt = [];
        T_skew = [];
        T_range = [];
        T_mean = [];
        T_median = [];
        T_var = [];
        T_std = [];
        T_cv = [];
        T_prc25 = [];
        T_prc75 = [];
        T_Trimean = [];
        T_shannon = [];
        T_entropy = [];
        T_Contrast= [];
        T_Correlation= [];
        T_Energy= [];
        T_Homogeneity= [];
        T_shuffle_index = [];
        T_G = [];
        T_GPP = [];
        


        %import disturebance reference cubes
        % mu = 0.01;
        % alpha = 1.0;
        % beta = 0.03;
        in_cubedir = '/Net/Groups/BGI/scratch/swang/Data/D_cube/intensity/';
        f_flag = sprintf('D_cube_Intenstiy_mu_%.3f_alpha_%.3f_beta_%.3f_on_*.mat',mu,alpha,beta);
        f = dir(fullfile(in_cubedir,f_flag));
        filepath = f.name;
        cubedir = fullfile(in_cubedir,filepath);
        
        load(cubedir)
        D(D>=1)=1;% all intensity values greater than 1 are 1!!!
        
        % now shuffle the D cubes.
        randIndex = randperm(size(D,3));
        D = D(:,:,randIndex);
        gpp_curves = zeros(N,11);
        gppFins = zeros(X,Y,11);
        BioFins = zeros(X,Y,11);

        for i = 1:11 % shuffle loop
            
            shuffle_index = i-1;
            %shuffle again
            randIndex = randperm(size(D,3));
            D_new = D(:,:,randIndex);

            %run Dynamic model
            [cVeg_D,gpp_D,npp_D,loss_D] = cdyn_DisT(c0,fGain,fLoss,a,b,e,k,N,D_new);
            
            %final 10-year steady biomass 
            Bio_10 = mean(cVeg_D(:,:,end-9:end),3);%average Biomass for last 10 year
            clear cVeg_D npp_D loss_D
            
            %gpp info
            gppFin = gpp_D(:,:,end-1);
            gpp_D = reshape(gpp_D,[],N);
            gpp_curve = mean(gpp_D,1);
            gpp = gpp_curve(end-1);

            %%statistics
            bio_unscale = filloutliers(Bio_10,'nearest','mean'); %very important, exclude the outliers with nearest non-outliers
            %---------------------------------------save No.1
            %save the biomass
            BioFins(:,:,i) = bio_unscale;
            
            % 11 statistics
            Bio = bio_unscale(:);
            kurt_bio = kurtosis(Bio);
            skew_bio = skewness(Bio);
            range_bio = prctile(Bio,90)-prctile(Bio,10);
            mean_bio = mean(Bio);
            median_bio = median(Bio);
            var_bio = var(Bio); %varience
            std_bio = std(Bio);
            cv_bio = 100*std_bio./mean_bio;
            prc25_bio = prctile(Bio,25);
            prc75_bio = prctile(Bio,75);
            clear Bio
            Trimean = 0.25*prc25_bio + 0.5*median_bio + 0.25*prc75_bio; 
            %entropy
            bioarray = uint8(rescale(bio_unscale, 0, 255));
            clear bio_unscale
            e = entropy_scale(bioarray);
            %shannon
            tb = tabulate(bioarray(:));  % number of values
            [s,VarS]=index_SaW(tb(:,2));
            %texture
            glcms = graycomatrix(bioarray);
            clear bioarray
            stats = graycoprops(glcms,{'contrast','Correlation','Energy','homogeneity'});
            clear glcms
            Contrast= stats.Contrast;
            Correlation= stats.Correlation;
            Energy= stats.Energy;
            Homogeneity= stats.Homogeneity;

            %---------------------------------------save No.3
            %save the xls
            %% add variables values to table
            T_mu(end+1) = mu;
            T_Is(end+1) = beta;
            T_alpha(end+1) = alpha;           
            T_kurt(end+1) = kurt_bio;
            T_skew(end+1) = skew_bio;
            T_range(end+1) = range_bio;
            T_mean(end+1) = mean_bio;
            T_median(end+1) = median_bio;
            T_var(end+1) = var_bio;
            T_std(end+1) = std_bio;
            T_cv(end+1) = cv_bio;
            T_prc25(end+1) = prc25_bio;
            T_prc75(end+1) = prc75_bio;
            T_Trimean(end+1) = Trimean;
            
            T_shannon(end+1) = s;
            T_entropy(end+1) = e;
            T_Contrast(end+1)= Contrast;
            T_Correlation(end+1)= Correlation;
            T_Energy(end+1)= Energy;
            T_Homogeneity(end+1)= Homogeneity;
            T_shuffle_index(end+1) = shuffle_index;
            T_G(end+1) = GPP_index;
            T_GPP(end+1) = gpp;

            %---------------------------------------save No.3
            %save the GPP
            gpp_curve = uint32(gpp_curve.*10);
            gppFin= uint32(gppFin.*10);
            gpp_curves(:,i) = gpp_curve;
            gppFins(:,:,i) = gppFin;
            
            
            clear D_new cVegFin Bio_10 gpp_curve gpp_10 gppFin cVegFin
            
        end

        %flag
        time2 = string(datetime('now'));
        lag = sprintf('mu_%.3f_a_%.3f_Beta_%.3f_GPP_%.2f_k_%.2f_are done at %s',mu,alpha,beta,GPP_index,Lb,time2);
        disp(flag)
        
        %save gpp info file
        %save(out_gppfullname,'gpp_curves','gppFins')
        %save BioFin info file 
        %save(out_cvegfullname,'BioFins')
        %save xls
        T = table(T_mu',T_Is',T_alpha',T_mean',T_median',T_range',T_var',T_std',T_cv',T_skew',T_kurt',T_prc25',T_prc75',T_Trimean',T_shannon',T_entropy',T_Contrast',T_Correlation',T_Energy',T_Homogeneity',T_shuffle_index',T_G',T_GPP');
        T.Properties.VariableNames = {'mu','beta','alpha','mean_bio','median_bio','range_bio','var_bio','std_bio','cv_bio','skew_bio','kurt_bio','prc25_bio','prc75_bio','Trimean','shannon','entropy','contrast','correlation','energy','homogeneity','Shuffle_index','G','GPP'};
        writetable(T,xlsfile)
    end
end






function [c,g,n,l] = cdyn_DisT(c0,fGain,fLoss,a,b,e,k,N,Pixel_index) % a b e k should be matrix with size of size(c0)
    c           = NaN([size(c0),N]);
    g           = c;
    n           = c;
    l           = c;
    c(:,:,1)      = c0;
    
    
    for i = 2:N
        Cveg_init     = c(:,:,i-1).*(1-Pixel_index(:,:,i-1));
        g(:,:,i-1)    = fGain(Cveg_init,a,b,e);
        n(:,:,i-1)    = g(:,:,i-1).*0.5;
        l(:,:,i-1)    = fLoss(Cveg_init,k);
        c(:,:,i)      = Cveg_init + n(:,:,i-1) - l(:,:,i-1); 
    end  

end

function [x] = entropy_scale(array)
    array = uint8(rescale(array, 0, 100));
    x = entropy(array);
end

% calculate Shannon-Index
function [H,VarH]=index_SaW(A,base)
    if nargin < 2
          base = 2; %default log2
    end
    [r,c] = size(A);
    warning('off','MATLAB:dispatcher:InexactMatch')
    if r==1
        A=A';
    end
    A(find(isnan(A))) = 0;
    n = sum(A); 
    N = repmat(n,r,1);
    B = A./N.*log(A./N)./log(base);
    B(find(isnan(B))) = 0;
    H = -sum(B);
    %variance of the Shannon-Wiener index
    aa = A./A;
    aa(find(isnan(aa))) = 0;
    S = sum(aa);
    VarH = (sum(N./A.*(B.^2))-(sum(B)).^2)./n+(S-1)./(2*n.^2);
end