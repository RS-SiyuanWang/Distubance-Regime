
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
    %domain_check = false(Domain_width,Domain_width);

   
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
                check_flag = max(domain(:));
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


