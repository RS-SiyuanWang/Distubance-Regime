
% DisturbanceRegimeParameters : mu, alpha, Iz(currently only slope);
% Biomass Statistics: mean, median, range, skewness, kurtosis.
function cube_intensity_batch(n)
    %% mu_n 1:9
    %% alpha 1:17

    load('parameter.mat')
    mu = parameter(n,2);
    alpha = parameter(n,3);

    I_slope = 0.078514; % parameter1 for Iz
    I_bias  = 0.22684; % parameter2 for Iz
    X = 1000; %window width of analysis
    Y = 1000; %window height of analysis
    arrsize = [X Y];
    

    % mu_v = [0.01:0.005:0.05];
    % alpha_v = [1.0:0.05:1.8];
    beta_v = [0.03:0.01:0.09 0.1:0.05:0.25 0.3:0.1:0.5];
    

    % diary 
    diary_folder = './Diary';
    if ~exist(diary_folder,'dir')
        mkdir(diary_folder)
    end
    diary_name = strcat('Diary_',datestr(now,'yyyy_mm-dd_HH:MM:SS'),'__Mu__',num2str(mu),'___Alpha__',num2str(alpha));
    diary(fullfile(diary_folder,diary_name))



    % load raw binary mat
    indir = '/Net/Groups/BGI/scratch/swang/Data/D_cube/raw/';    
    filename = sprintf('D_Cube_mu_%.3f_alpha_%.3f_20220419_1000x1000.mat',mu,alpha);
    load(fullfile(indir,filename))
    
    %slopes = [0.03:0.01:0.09 0.1:.1:.5];
   
    
    for x = 1:numel(beta_v)
        I_slope = beta_v(x);
        
        check = [1,2,3,10,100,200];
        for i=1:size(D_cube,3)
            D_i = D_cube(:,:,i);
            TotalPixel = sum(D_i(:));
            [I_array,D_map,num] = IntensityForAreas(I_slope,I_bias,D_i);
            D(:,:,i) = I_array;
            time = string(datetime('now'));
            text = sprintf('NO.%d_____mu_%.3f_a_%.3f at %s, disturbed pixels are %d',i,mu,alpha,time,TotalPixel);
            if ismember(i,check) 
                disp(text);
            end
        end

        % save outfile 
        outname = sprintf('D_cube_Intenstiy_mu_%.3f_alpha_%.3f_beta_%.3f_on_20220421_1kx1k_Disturbed_%d.mat',mu,alpha,I_slope,TotalPixel);
        outdir = '/Net/Groups/BGI/scratch/swang/Data/D_cube/intensity/';
        outdir = fullfile(outdir,outname);
        save(outdir,'D')

        time2 = string(datetime('now'));
        flag = sprintf('%s is done at %s',outdir,time2);
        disp(flag)
        
    end
    diary off
end

function [I_array,D_map,num] = IntensityForAreas(a,b,Disturbance_Binary_map)
    
    D_map = Disturbance_Binary_map;
    [L,num] = bwlabel(Disturbance_Binary_map,4);
    areas = regionprops(L,'Area');
    %status = regionprops(L,'BoundingBox');
    %centroid = regionprops(L,'Centroid');
    Iz = @(a,b,area) a.*log10(100.*area)+b;
    I_array = L;
    for i=1:num
        % rectangle('position',status(i).BoundingBox,'edgecolor','r');
        area = areas(i,1).Area(1);
        I = Iz(a,b,area);
        I_array(L == i) = I;
        %text(centroid(i,1).Centroid(1,1),centroid(i,1).Centroid(1,2), num2str(area),'Color', 'r') 
    end
end
