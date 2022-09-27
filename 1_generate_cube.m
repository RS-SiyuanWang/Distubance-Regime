

% % test the result
% [D_cube] = generate_cube_(1);
% domain = D_cube(:,:,1);
% figure(1);
% imagesc(domain)
% colormap('jet')
% colorbar

function generate_cube(n)

    %tic
    load('parameter.mat')
    mu = parameter(n,2);
    a = parameter(n,3);
    %beta = index(n,4);
    
    DomainWidth = 1000;
    DomainHeight = 1000;
    CubeDepth = 200;
    Domains = DomainWidth .* DomainHeight; % Number of Domain
 
    % start to diary
    diary_folder = './Diary';
    if ~exist(diary_folder,'dir')
        mkdir(diary_folder)
    end
    diary_name = strcat('Diary_',datestr(now,'yyyy_mm-dd_HH'),'_index_',num2str(n));
    diary(fullfile(diary_folder,diary_name))

    tic
    error = 1;
    while error > 0.001
        [TotalPixel, GapSizes, GapNums] = Generate_GapNums(DomainWidth,mu,a);
        error = abs(TotalPixel-Domains.*mu)./Domains.*100;  % error means the pencentage of difference between expected mu and generated mu
    end
        
    D_cube = zeros(DomainHeight, DomainWidth, CubeDepth);

    % xls file for the disturbed numbers
    T_mu=[];
    T_a =[];
    T_theoretic = [];
    T_prescribed = [];
    T_generated = [];
    T_depth = [];
        
    %output the checkpoints at 3 snapshots
    check = [1,10,100,200];
    for i = 1:CubeDepth
        %t1=cputime;
        [domain, ExistNum, XPlist, YPlist, Wlist, Llist, DisPixels] = randomly_place_rect(GapNums,GapSizes,DomainWidth);
        D_cube(:,:,i) = domain;
        if ismember(i,check)
            time = string(datetime('now'));
            text = sprintf('mu_%.3f_a_%.3f_at %s, disturbed pixels are %d/%d',mu,a,time,TotalPixel,DisPixels);
            disp(text);
        end

        T_mu(end+1)=mu;
        T_a(end+1) =a;
        T_theoretic(end+1) = Domains.*mu;
        T_prescribed(end+1) = TotalPixel;
        T_generated(end+1) = DisPixels;
        T_depth(end+1) = i;
    end
        
    toc
    

    outdir = '/Net/Groups/BGI/scratch/swang/Data/D_cube/raw/';
    outname = sprintf('D_Cube_mu_%.3f_alpha_%.3f_20220419_1000x1000.mat',mu,a);
    save(fullfile(outdir,outname),'D_cube')
    disp(outname);

    xlsdir = '/Net/Groups/BGI/scratch/swang/Data/D_cube/raw/xls/';
    xlsname = sprintf('D_Cube_mu_%.3f_alpha_%.3f_20220419_1000x1000.xls',mu,a);
    xlsfile = fullfile(xlsdir,xlsname);
        
    T = table(T_mu',T_a',T_theoretic',T_prescribed',T_generated',T_depth');
    T.Properties.VariableNames = {'mu','alpha','theoretic_area','prescribed_area','generated_area','depth'};
    writetable(T,xlsfile)

    diary off
    
    
end

function [TotalPixel, GapSizes, GapNums] = Generate_GapNums(DomainWidth,mu,a)
    % (x) GapSize needs to be modified for different sizes of domain
    GenerateGapNums = @(A,a,X) A*X.^(-a);
    TotalPixel = 0;
    GapNums =[];
    GapSizes = [];
    sum_ = 0; % initialization of sum to generate proportionality
    Max_exp = floor(log2(DomainWidth.*DomainWidth)); % maximum exponent of 2
    for i = 0:Max_exp
        GapSize = 2^i;
        GapSizes = [GapSizes GapSize];
    end
    
    for i = 1: numel(GapSizes)
        Zi = GapSizes(i);
        sum_ = sum_ + Zi .* Zi.^(-a) ;
    end
    D = DomainWidth .* DomainWidth;
    A = D .* mu ./ sum_; % proportionality parameter A 
    GapNums = GenerateGapNums(A,a,GapSizes);
    %GapNums = round(GapNums);
    % pesudo-random process to GapNums, make it to integer value.
    for i = 1:numel(GapNums)
        Gap_i = GapNums(i);
        remainder = rem(Gap_i,1);
        if rand(1) < remainder
            GapNums(i) = ceil(Gap_i);
        else
            GapNums(i) = fix(Gap_i);
        end
    end
    for i =1:numel(GapNums)
        TotalPixel = TotalPixel + GapNums(i).*GapSizes(i);
        i = i+1;
    end
end



function [domain, ExistNum, XPlist, YPlist, Wlist, Llist, TotalPixels] = randomly_place_rect(GapNums, GapSizes, Domain_width)

    squarePixels = @(Width) zeros(Width, Width) + 1;
    rectPixels = @(Width, Length) zeros(Length, Width) + 1;
    Gaps_N = sum(GapNums(:)); % the amounts of total gaps
    % Gaps_N = sum(GapNums, 'all'); % the amounts of total gaps
    Width_vector = zeros(1, Gaps_N); % initialization of width vector contains every gaps
    Length_vector = zeros(1, Gaps_N);
    Width_vector = Width_vector(:);
    Length_vector = Length_vector(:);

    Widths = [];
    Lengths = [];

    % intensity anonymous function, a for beta, b for the intercept, equals to 0.22684
    % beta_intercept = 0.22684;
    % Iz = @(a,b,area) a.*log10(100.*area)+b;

    for i = 1:numel(GapSizes)

        if rem(i, 2) == 0
            % Serial Number is Even, the corresponding square root is not an
            % integer;
            width = sqrt(GapSizes(i - 1));
            length = GapSizes(i) ./ width;
            Widths(i) = width;
            Lengths(i) = length;
        else
            width = sqrt(GapSizes(i));
            Widths(i) = width;
            Lengths(i) = width;
        end

    end

    % generate vectors of widths and lengths to plot
    N = 1;

    for i = 1:numel(GapNums)
        GapNums_i = GapNums(i);
        Width_vector(N:N + GapNums_i - 1) = Widths(i);
        Length_vector(N:N + GapNums_i - 1) = Lengths(i);
        N = N + GapNums_i;
    end

    % define the domain matrix
    Xmax = Domain_width;
    Ymax = Domain_width;

    
    Widths_v = Width_vector;
    Lengths_v = Length_vector;
    result = Widths_v .* Lengths_v;
    expect = GapNums .* GapSizes;

    if sum(result(:)) - sum(expect(:)) == 0
    %if sum(Widths_v .* Lengths_v, 'all') - sum(GapNums .* GapSizes, 'all') == 0
        %text = ['Total Pixel generated succeed, ', num2str(sum(Widths_v .* Lengths_v, 'all')), ' in total'];
        text = ['Total Pixel generated succeed, ', num2str(sum(expect(:))), ' in total'];
        %disp(text)
    else
        disp('Gap pixels number is wrong')
    end

    % Initialization Parameters
    ExistNum = 0; % the numbers of the existing squares
    XPlist = []; % the x position of the center of existing squares
    YPlist = []; % the y position of the center of existing squares
    Wlist = []; % the widthes of existing squares
    Llist = []; % the lengthes of existing squares

    domain = false(Domain_width,Domain_width);
    domain_check = false(Domain_width,Domain_width);

   
    %place each event
    
    for i=1:Gaps_N
        
        w = Widths_v(end - i + 1); % from large to small
        l = Lengths_v(end - i + 1);

        check_flag = 5;
        while check_flag > 1
            xp = randi([1,Xmax]);
            yp = randi([1,Ymax]);
            %disp(xp)
            if xp + w ./ 2 <= Xmax && xp - w ./ 2 >= 0 && yp + l ./ 2 <= Ymax && yp - l ./ 2 >= 0 % no edge connecting
                I0 = false(Domain_width,Domain_width);
                I0(xp,yp) = 1;
                SE = strel(rectPixels(l, w));
                Event = imdilate(I0, SE);
                domain_check = domain + Event;
                check_flag = max(domain_check(:));
            else
                continue
            end 
        end 

        ExistNum = ExistNum + 1;
        Wlist(ExistNum) = w;
        Llist(ExistNum) = l;
        XPlist(ExistNum) = xp;
        YPlist(ExistNum) = yp;

        % e_size = w.*l;
        % *intensity = Iz(beta,beta_intercept,e_size);
        % Event_intensity = Event .* intensity;
        domain = domain + Event;
    end

    prescribed = sum(expect(:));
    TotalPixels = sum(domain(:));
    flag = sprintf('Total %8d events generated, %d/%d generated/prescibed',ExistNum,TotalPixels,prescribed);
    %disp(flag)
end


