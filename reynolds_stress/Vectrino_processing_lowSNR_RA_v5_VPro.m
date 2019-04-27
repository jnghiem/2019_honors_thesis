%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vectrino_processing_lowSNR
%The input is the .dat file created by converting the .vno file in the
%Vectrino software.

%Modified from Vectrino_processing_final by Laurel on 8/10/2011. 
%Modified to use Vectrino Profiler Data by Rachel Allen (rachelallen@berkeley.edu) May 2014.
% Updated with despiking_on and cell_select options.  Rachel Allen, Sept 2014

%This code does not use an SNR filter but instead eliminates most bad datapoints
%through the threshold despiking algorithm. It is not recommended to use
%this code when interested in turbulent statistics, as the noise in the
%data may be confounded with actual turbulence. However, the data can be
%used to calculate means. 

%This code also differs from Vectino_processing_final in that it does not
%replace bad datapoints at all (since there are many more bad datapoints
%than good datapoints). Derivative and second derivative surrogates in the threshold
%despiking algortihm are calculated using actual times but are then
%divided by 600 to meet the criteria (according to Goring and Nikora)
%that the standard deviation of velocity be approximately the same order of
%magnitude as the standard deviation of the second derivative surrogate.
%(In contrast, in Vectrino_processing_final, derivative surrogates are not
%constructed according to time. That code only works when all points are
%equally spaced.) What this means is that when applying
%Vectrino_processing_lowSNR to a new dataset with very different velocities
%from the Everglades, the user should check the magnitudes of these
%standard deviations and adjust the constant multiplier if needed before
%implementing the threshold despiking algorithm.

%Note that the code for using a correlation filter and low-SNR filter is
%still present and may be turned on or off by commenting or uncommenting
%those lines.

%NOTE: To continue to avoid replacing bad datapoints, as recommended, make
%sure afterGN is set to 'y' below!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('C:\Users\rneuhausler\Everglades Research')
load('filename_VPro')


for f = 1:length(filename_VPro)
    
clearvars -except f*
close all

    
    f1 = filename_VPro{f};
    filename = [f1(1:end-4),'.mat'];

cd('C:\Users\rneuhausler\Everglades Research\Raw_data_2013_November_Vectrino-II')

%USER INPUT SECTION
%filename = 'C:\Users\Rachel\Projects\Berkeley\Research\VectrinoII\2013-11_Everglades\01_raw\Post5-4-avg20131106101043.dat'; %This is the input .dat file from the Vectrino software.
%filename = 'C:\Users\Rachel\Projects\Berkeley\Research\VectrinoII\2013-11_Everglades\01_raw\post5-2-avg20131105123253.dat'; %This is the input .dat file from the Vectrino software.
%filename = 'C:\Users\Rachel\Projects\Berkeley\Research\VectrinoII\2013-11_Everglades\01_raw\RS1uPreFlow.307.13.Vectrino Profiler1.00001.mat'; % Vectrino Profiler Data
%filename = 'C:\Users\Rachel\Projects\Berkeley\Research\VectrinoII\2013-11_Everglades\01_raw\post2-5-.311.13.Vectrino Profiler.00002.mat'; % Vectrino Profiler Data
% data Rosanna selected as "good"
%filename = 'C:\Users\Rachel\Projects\Berkeley\Research\VectrinoII\2013-11_Everglades\01_raw\Post2-4-avg20131106115505.dat';



%filename = 'post1-5-.311.13.Vectrino Profiler.00000.mat';
%filename = 'post1-5-.311.13.Vectrino Profiler.00001.mat';
%filename = 'post1-5-.311.13.Vectrino Profiler.00002.mat';
%filename = 'post1-5-.311.13.Vectrino Profiler.00003.mat';
%filename = 'post1-5-.311.13.Vectrino Profiler.00004.mat';
%filename = 'post1-5-.311.13.Vectrino Profiler.00005.mat';

%filename = 'post2-5-.311.13.Vectrino Profiler.00000.mat';
%filename = 'post2-5-.311.13.Vectrino Profiler.00001.mat';
%filename = 'post2-5-.311.13.Vectrino Profiler.00002.mat';
%filename = 'post2-5-.311.13.Vectrino Profiler.00003.mat';

%filename = 'post3-5-.311.13.Vectrino Profiler.00000.mat';
%filename = 'post3-5-.311.13.Vectrino Profiler.00001.mat';
%filename = 'post3-5-.311.13.Vectrino Profiler.00002.mat';
%filename = 'post3-5-.311.13.Vectrino Profiler.00003.mat';

%filename = 'post4-5-.311.13.Vectrino Profiler.00000.mat';
%filename = 'post4-5-.311.13.Vectrino Profiler.00001.mat';
%filename = 'post4-5-.311.13.Vectrino Profiler.00002.mat';






%to_save = 'C:\Users\Rachel\Projects\Berkeley\Research\VectrinoII\2013-11_Everglades\02_processed\Post5-4-avg20131106101043.mat'; %File to save output
%to_save = 'C:\Users\Rachel\Projects\Berkeley\Research\VectrinoII\2013-11_Everglades\02_processed\post1-5-.311.13.Vectrino Profiler.00002.mat'; %File to save output
%to_save = 'C:\Users\Rachel\Projects\Berkeley\Research\VectrinoII\2013-11_Everglades\02_processed\Post2-4-avg20131106115505.mat';

to_save = filename;

n_columns = 19; %This is the number of columns in the .dat file. It will vary depending on the options selected in the conversion process.
time_included = 'y'; %'y' if the data file includes times, 'n' if not. (This will either be column 1 or 2. Contrary to what the hdr file says, time is in seconds!)
afterGN = 'n'; %'y' if bad data points are replaced through cubic interpolation after the Goring + Nikora procedure has converged; 'n' if replacement occurs within the procedure. 'n' is default, but change to 'y' if the replacements are diverging.
corr_cutoff = 40; %Percent correlation cutoff for good data
SNR_cutoff = 5; %SNR (dB) cutoff for good data
pct_diff_cutoff = 50; %Datapoints with v_z1 being more than this percent different from v_z2 will be excluded, unless they are within the absolute difference cutoff below. Default is 50.
abs_diff_cutoff = 0.002; %m/s. z-direction datapoints that exceed the percent difference cutoff above but are within this absolute difference will be retained. Default is 0.002.
threshold_pd_cutoff = 0.1; %When the size of the ellipse used to exclude points in the threshold despiking algorithm changes by less than this amount, convergence is achieved. Default is 0.1. 
calc_turbulence = 'y'; %Set this to yes when you want to calculate and save turbulence statistics.
interp = 1; %Interpolate bad values?  1 = yes, 0 = no.
isProfiler = 1; % Is the data from a Vectrino Profiler?  1 = yes, 0 = no.
despiking_on = 0; % Use Goring-Nikora despiking algorithm? 1 = on, 0 = off
cell_select = NaN; % select a specific cell (1 through nCells) in Vectrino profiler data.  NaN = use all cells

if isProfiler % for Profiler Data
% load VII data
load(filename)
%select out the data that we want
v_vII(:,:,1) = Data.Profiles_VelX;
v_vII(:,:,2) = Data.Profiles_VelY;
v_vII(:,:,3) = Data.Profiles_VelZ1;
v_vII(:,:,4) = Data.Profiles_VelZ2;
SNR_vII(:,:,1) = Data.Profiles_SNRBeam1;
SNR_vII(:,:,2) = Data.Profiles_SNRBeam2;
SNR_vII(:,:,3) = Data.Profiles_SNRBeam3;
SNR_vII(:,:,4) = Data.Profiles_SNRBeam4;
corr_vII(:,:,1) = Data.Profiles_CorBeam1;
corr_vII(:,:,2) = Data.Profiles_CorBeam2;
corr_vII(:,:,3) = Data.Profiles_CorBeam3;
corr_vII(:,:,4) = Data.Profiles_CorBeam4;
t_vII = Data.Profiles_TimeStamp;
nCells = Config.nCells; % number of cells

% CREATE CELL ARRAYS TO STORE THE CLEANED RESULTS
v_clean = cell(nCells,1);
t_clean = cell(nCells,1);

% record the total percent of data removed
pct_rm_clean = cell(nCells,1);

else % for Point ADV data
% %READ IN DATA FILE
fid = fopen(filename);
nCells = 1;
switch n_columns
    case 19
        if time_included == 'y'
            A = textscan(fid, '%n%u32%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n');
            t =transpose(0:0.1:(length(A{1})-1)*0.1)/60; %time, min
        elseif time_included == 'n'
            A = textscan(fid, '%u32%u32%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n');
            t = [];
        else
            error('Enter "y" or "n" for "time_included" and rerun.')
        end
        counter = A{2}; %Number of the datapoint
        v = [A{4}, A{5}, A{6}, A{7}]; %vx, vy, vz, vz2
        SNR = [A{12}, A{13}, A{14}, A{15}]; %SNR, dB, beams 1-4
        corr = [A{16}, A{17}, A{18}, A{19}]; %correlation, %, beams 1-4
    case 20
        A = textscan(fid, '%u32%n%u32%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n');
        t = A{2}/60; %time, min
        counter = A{3}; %Number of the datapoint
        v = [A{5}, A{6}, A{7}, A{8}]; %vx, vy, vz, vz2
        SNR = [A{13}, A{14}, A{15}, A{16}]; %SNR, dB, beams 1-4
        corr = [A{17}, A{18}, A{19}, A{20}]; %correlation, %, beams 1-4
    case 18
        A = textscan(fid, '%u32%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n');
        t = [];
        counter = A{1}; %Number of the datapoint
        v = [A{3}, A{4}, A{5}, A{6}]; %vx, vy, vz, vz2
        SNR = [A{11}, A{12}, A{13}, A{14}]; %SNR, dB, beams 1-4
        corr = [A{15}, A{16}, A{17}, A{18}]; %correlation, %, beams 1-4
end
fclose(fid);
clear A

pct_rm = NaN; % store the percent of data removed

end

%DETERMINE WHETHER THE DATA RECORD NEEDS TO BE TRUNCATED
%figure, plot(t,v)
% response = 'bad';
% while strcmp(response, 'bad')
%     truncate = input('Truncate data record? \n', 's');
%     switch truncate
%         case {'y', 'yes', 'Y', 'Yes'}
%             startt = input('\nEnter start time for good data or push enter to keep current start time\n');
%             endt = input('\nEnter end time for good data or push enter to keep current start time\n');
%             if isempty(startt), startt = min(t); end
%             if isempty(endt), endt = max(t); end
%             if ~isempty(t) 
%                 v = v(t<=endt, :); corr = corr(t<=endt); SNR = SNR(t<=endt); t = t(t<=endt); 
%                 v = v(t>=startt,:); corr = corr(t>=startt); SNR = SNR(t>=startt); t = t(t>=startt);
%             else
%                 v = v(counter<=endt, :); corr = corr(counter<=endt); SNR = SNR(counter<=endt); counter = counter(counter<=endt);
%                 v = v(counter>=startt,:); corr = corr(counter>=startt); SNR = SNR(counter>=startt); counter = counter(counter>=startt);
%             end
%             response = 'good';
%         case {'n', 'no', 'N', 'No'}
%             response = 'good';
%         otherwise
%             response = 'bad';
%             disp('Enter yes or no')
%     end
% end

% select the cells to display
if isnan(cell_select)
    c_start = 1;
    c_end = nCells;
else
    c_start = cell_select;
    c_end = cell_select;
end
%%

for c = c_start:c_end % treat each depth cell independently
    
    if isProfiler
    v = squeeze(v_vII(:,c,:));
    SNR = squeeze(SNR_vII(:,c,:));
    corr = squeeze(corr_vII(:,c,:));
    t = t_vII; % refresh timestamp
    end
    tmin = min(t); 
    tmax = max(t); % store limits for plotting
    
    % find length of data at start
    length_v_start = length(v);
    
    figure(6), 
    ha(1) = subplot(3,1,1);
    plot(t,v,'.-'), title(['cell = ',num2str(c),', before filtering'])%,pause
    xlim([tmin,tmax])
    
exclude = [];

%REMOVE DATA THAT DOESN'T MEET THE PERCENT CORRELATION CUTOFF
[I, J] = find(corr<corr_cutoff);
exclude = unique(I); %Indices of data array in which at least one of the beams didn't meet the percent correlation cutoff.
clear I J

% % REMOVE DATA THAT DOESN'T MEET THE SNR CUTOFF
% [I,J] = find(SNR<SNR_cutoff);
% exclude = unique([exclude; I]); %Indices of the data array in which at least one of the beams didn't meet the percent correlation OR SNR cutoff.
% clear I J
% 
%REDUNDANT Z FILTER: REMOVE THE DATA WITH LARGE DISCREPANCIES BETWEEN V_Z1 AND V_Z2
pct_diff = abs((v(:,3)-v(:,4))./((v(:,3)+v(:,4))/2))*100; %Percent difference between the two z-direction velocities
abs_diff = abs(v(:,3)-v(:,4)); %Absolute difference between the two z-direction velocities
exclude = unique([exclude; intersect(find(pct_diff>pct_diff_cutoff), find(abs_diff>abs_diff_cutoff))]);
clear pct_diff abs_diff

%REMOVE BAD DATA IDENTIFIED ABOVE
include = setdiff(1:size(v,1), exclude); %Indices of good data
v = v(include,:); t = t(include); %exclude bad datapoints
exclude = []; %Reset
oldv = v; %Within the loop below, v will be modified by subtracting the mean. This saves the original v. Only for replacement after G+N
oldt = t;
if afterGN == 'y'
    oldinclude = 1:size(v,1); %Only for replacement after G+N
else
    Exclude = [];
end
addback = [0 0 0 0]; %what needs to be added back to v at the end of the threshold despiking procedure
clear include
figure(6)
ha(2) = subplot(3,1,2);
plot(t,v,'.-'), title(['cell = ',num2str(c),', after SNR, Corr, and Z filters'])%,pause
    xlim([tmin,tmax])
figure(4), close figure 4

%PERFORM THRESHOLD DESPIKING ALGORITHM (GORING AND NIKORA)
if despiking_on
pctdiff = 1; olda = [0 0 0 0]; oldb = [0 0 0 0]; %initialize
niter = 1;
size_exclude = 5;
while pctdiff>threshold_pd_cutoff && niter < 100 && size_exclude > 4
      if size(v,1) < 10
        msg{c} = ['too many data points removed, check data quality and removal parameters and re-run']
        break
      end
      
    %Subtract the mean
    mean_row = mean(v, 1);
    v = v-repmat(mean_row, size(v,1),1);
    addback = addback+mean_row;
    
    %Calculate derivative surrogates:
    del_v = (v(3:size(v,1), :)-v(1:size(v,1)-2, :))./repmat(t(3:size(v,1))-t(1:size(v,1)-2), 1, 4)/600; %Profile at each end is eliminated.
    t_v = t(2:size(v,1)-1);
    del2_v = (del_v(3:size(del_v,1),:)-del_v(1:size(del_v,1)-2,:))./repmat(t_v(3:size(del_v,1))-t_v(1:size(del_v,1)-2),1,4)/600;
    lambda_v = std(v)*sqrt(2*log(size(v,1)));
    lambda_del_v = std(del_v)*sqrt(2*log(size(del_v,1)));
    lambda_del2_v = std(del2_v)*sqrt(2*log(size(del2_v,1)));
    theta = atan(sum(v(3:size(v,1)-2,:).*del2_v)./sum(v(3:size(v,1),:).^2));
    b = ((lambda_del2_v.^2-lambda_v.^2.*(tan(theta)).^2)./((cos(theta)).^2-(tan(theta)).^2.*(sin(theta)).^2)).^0.5;
    a = (lambda_v.^2/(cos(theta)).^2-b.^2.*(tan(theta)).^2).^0.5;
    th = 0:pi/180:2*pi;

    %Eliminate points that are outside the ellipses in each projection
    tags34 = ones(size(v,1),2); %This matrix will be used to determine which of the two z-direction velocities (or an average of the two) should be used
    tags12 = ones(size(v,1),1); 
    for beam = 1:4
        [data_theta_1 data_r_1] = cart2pol(v(2:size(v,1)-1, beam), del_v(:,beam));
        [data_theta_2 data_r_2] = cart2pol(del_v(2:size(del_v,1)-1, beam), del2_v(:,beam));
        [data_theta_3 data_r_3] = cart2pol(v(3:size(v,1)-2, beam), del2_v(:,beam));
        r1_exp = lambda_v(beam)*lambda_del_v(beam).*((lambda_del_v(beam).*cos(data_theta_1)).^2+(lambda_v(beam).*sin(data_theta_1)).^2).^-0.5;
        r2_exp = lambda_del_v(beam)*lambda_del2_v(beam).*((lambda_del2_v(beam).*cos(data_theta_2)).^2+(lambda_del_v(beam).*sin(data_theta_2)).^2).^-0.5;
        r3_exp = (sqrt(2)*a(beam)*b(beam)*((b(beam)^2-a(beam)^2).*cos(2*data_theta_3-2*theta(beam))+a(beam)^2+b(beam)^2).^0.5)./((b(beam)^2-a(beam)^2).*cos(2*data_theta_3-2*theta(beam))+a(beam)^2+b(beam)^2);
        to_exclude1 = find(data_r_1>r1_exp)+1; to_keep1 = setdiff(2:size(v,1)-1, to_exclude1);
        to_exclude2 = find(data_r_2>r2_exp)+2; to_keep2 = setdiff(3:size(v,1)-2, to_exclude2);
        to_exclude3 = find(data_r_3>r3_exp)+2; to_keep3 = setdiff(3:size(v,1)-2, to_exclude3);
        x1 = lambda_v(beam)*cos(th)*cos(0)-lambda_del_v(beam)*sin(th)*sin(0);
        y1 = lambda_v(beam)*cos(th)*sin(0)+lambda_del_v(beam)*sin(th)*cos(0);
        x2 = lambda_del_v(beam)*cos(th)*cos(0)-lambda_del2_v(beam)*sin(th)*sin(0);
        y2 = lambda_del_v(beam)*cos(th)*sin(0)+lambda_del2_v(beam)*sin(th)*cos(0);
        x3 = a(beam)*cos(th)*cos(theta(beam))-b(beam)*sin(th)*sin(theta(beam));
        y3 = a(beam)*cos(th)*sin(theta(beam))+b(beam)*sin(th)*cos(theta(beam));
        if beam == 2
        figure(1), clf
        subplot(2,2,1), plot(v(to_keep1, beam), del_v(to_keep1-1, beam), 'k.'), hold on, plot(v(to_exclude1, beam), del_v(to_exclude1-1, beam), 'g.'), plot(x1, y1, 'r-')
        subplot(2,2,2), plot(del2_v(to_keep2-2, beam), del_v(to_keep2-1, beam), 'k.'), hold on, plot(del2_v(to_exclude2-2, beam), del_v(to_exclude2-1, beam), 'g.'), plot(y2, x2, 'r-')
        subplot(2,2,3), plot(v(to_keep3, beam), del2_v(to_keep3-2, beam), 'k.'), hold on, plot(v(to_exclude3, beam), del2_v(to_exclude3-2, beam), 'g.'), plot(x3,y3, 'r-')
        title(sprintf('%s%d', 'Beam ', beam)), end
        switch beam
            case {1, 2}
                exclude = union(exclude, union(union(to_exclude3, union(to_exclude1, to_exclude2)), [1; 2; size(v,1); size(v,1)-1])); 
                tags12(exclude) = 0;
            case 3
                try
                    to_exclude = [union(to_exclude3, union(to_exclude1, to_exclude2)); 1; 2; size(v,1); size(v,1)-1];
                catch
                    to_exclude = [transpose(union(to_exclude3, union(to_exclude1, to_exclude2))); 1; 2; size(v,1); size(v,1)-1];
                end
                tags34(to_exclude, ones(size(to_exclude))) = 0; 
            case 4
                try
                    to_exclude = intersect(to_exclude, [union(to_exclude3, union(to_exclude1, to_exclude2)); 1; 2; size(v,1); size(v,1)-1]); %only exclude data if both v_z1 and v_z2 are bad.
                    tags34([union(to_exclude3, union(to_exclude1, to_exclude2)); 1; 2; size(v,1); size(v,1)-1], 2*ones(size(to_exclude))) = 0;
                catch
                    to_exclude = intersect(to_exclude, [transpose(union(to_exclude3, union(to_exclude1, to_exclude2))); 1; 2; size(v,1); size(v,1)-1]); %only exclude data if both v_z1 and v_z2 are bad.
                    tags34([transpose(union(to_exclude3, union(to_exclude1, to_exclude2))); 1; 2; size(v,1); size(v,1)-1], 2*ones(size(to_exclude))) = 0;
                end
                to_change = find(sum(tags34,2)==1 & tags12 == 1); %These are indices of when just one v_z datapoint is bad.
                if ~isempty(to_change)
                    v(to_change, 3:4) = [sum(tags34(to_change,:).*v(to_change, 3:4), 2), sum(tags34(to_change,:).*v(to_change, 3:4), 2)]; %Set bad v_z data equal to good v_z data at the same time
                end
                exclude = union(exclude, to_exclude); 
        end
    end
    clear del_v del2_v lambda_v lambda_del_v lambda_del2_v theta th beam data_theta_1 data_r_1 data_theta_2 data_r_2 data_thetha_3 data_r_3 r1_exp r2_exp r3_exp to_exclude1 to_exclude2 to_exclude3 x1 y1 x2 y2 x3 y3 to_exclude

    %REPLACE BAD DATA VALUES
    
%     %Replace remaining bad data with cubically interpolated values
    include = setdiff(1:size(v,1), exclude); %Indices of good data
    if afterGN == 'y'
        v = v(include, :); %Only for replacement after G+N
        oldinclude = oldinclude(include); %Only for replacement after G+N
        if ~isempty(t)
            t = t(include); %Truncate t, like v, to avoid data replacement errors at the ends of the data record
        else
            counter = counter(include);
        end
        fprintf('%s%d', 'Percent data exclusion is ', 100-length(oldinclude)/size(oldv,1)*100) %Only for replacement after G+N
    else
        v(exclude,:) = interp1(include, v(include, :), exclude, 'pchip', NaN);     
        [rownotnans colnotnans] = ind2sub(size(v), find(~isnan(v))); %Remove data values that would need to be extrapolated
        v = v(unique(rownotnans), :);
        if ~isempty(t)
            t = t(unique(rownotnans)); %Truncate t, like v, to avoid data replacement errors at the ends of the data record
        else
            counter = counter(unique(rownotnans));
        end
        clear colnotnans rownotnans
        Exclude = union(Exclude, exclude);
        fprintf('%s%d', 'Percent data exclusion is ', length(Exclude)/size(v,1)*100) %Only for replacement before G+N
    end
    size_exclude = length(exclude)
    exclude = []; %Reset
    clear include
    niter = niter+1
    clear tags
    pctdiff = max(abs(([a, b]-[olda, oldb])./[a, b]))*100
    olda = a; oldb = b;    
    figure(2), plot(t,v),title(['cell = ',num2str(c),', despiking process'])%,pause%, xlim([100-niter 2100-niter])
    figure(4), hold on, plot(niter, pctdiff, 'k.'), title(['cell = ',num2str(c)])
%    pause
end
clear a b olda oldb pctdiff


% %Add back the values that were subtracted from the v array in the threshold
% %despiking procedure. Only for replacement before G+K
if afterGN == 'n'
    v = v+repmat(addback, size(v,1),1);
else
    %Replace remaining bad data with cubically interpolated values. This section only for replacement after G+N
    v= oldv;
    exclude = setdiff(1:size(v,1), oldinclude); %Indices of good data
    v(exclude,:) = interp1(oldinclude, v(oldinclude, :), exclude, 'pchip', NaN);
    [rownotnans colnotnans] = ind2sub(size(v), find(~isnan(v))); %Remove data values that would need to be extrapolated
    v = v(unique(rownotnans), :);
    if ~isempty(t)
        t = t(unique(rownotnans)); %Truncate t, like v, to avoid data replacement errors at the ends of the data record
    else
        counter = counter(unique(rownotnans));
    end
    clear colnotnans rownotnans
end
end % end of DESPIKING

pct_rm = 1 - (length(v)/length_v_start); % find fraction of data removed, across all steps

if exist('msg') && size(msg{c},2) > 1
    v = NaN;
    t = NaN;
end

figure(6)
ha(3) = subplot(3,1,3);
plot(t,v,'.-'),title(['cell = ',num2str(c),', after all filters'])%,pause
    xlim([tmin,tmax])
    linkaxes(ha,'x')


if isProfiler
% store results
v_clean{c} = v;
t_clean{c} = t;
pct_rm_clean{c} = pct_rm;
clear v t SNR corr pct_rm tmin tmax
end

end

%%
% %COME UP WITH A SINGLE Z-DIRECTION VELOCITY THAT IS AN AVERAGE OF THE TWO 
% v = [v(:,1:2), mean(v(:,3:4), 2)]; %Replace the two v_z columns with the single one computed above

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %The section between these lines is specific to SharkTREX only
% 
% %ROTATE COLUMNS SO THAT X,Y, AND Z ARE ACTUALLY STREAMWISE (A.K.A.
% %BOAT-WISE), CROSS-STREAM, AND VERTICAL COMPONENTS WITH RESPECT TO THE
% %RIVER
% 
% v = [v(:,3), v(:,1), v(:,2)]; %Column order is now actual x, actual y, actual z
% 
% figure(3), plot(v)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATE TURBULENCE INTENSITIES AND REYNOLDS STRESS FOR THE WHOLE DATA
%RECORD IN THE FILE:
if calc_turbulence == 'y'
    if isProfiler

    vm = zeros(nCells,4);
    vp_clean = cell(nCells,1);
    vps_clean = vp_clean;
    vpvp_mean = zeros(4,4,nCells);
    
    for c = 1:nCells % treat each depth cell independently
             if  exist('msg') && size(msg{c},2) > 1
                break
             end
        vm(c,:) = mean(v_clean{c});
        vp_clean{c} = v_clean{c} - repmat(vm(c,:),size(v_clean{c},1),1);
        vps_clean{c} = vp_clean{c}.^2;
        for ii = 1:4
            for jj = 1:4
                vpvp_mean(ii,jj,c) = mean(vp_clean{c}(:,ii).*vp_clean{c}(:,jj),1);
            end
        end
    end
    else % point ADV
        if exist('msg') 
            vpvp_mean = NaN;
        else
x = mean(v(:,1));
y = mean(v(:,2));
z = mean(v(:,3));
z2 = mean(v(:,4));
vp = v-repmat([x,y,z, z2], size(v,1), 1);
vps = mean(vp.^2);
vpvp_mean = NaN(4);
for ii = 1:4
    for jj = 1:4
        vpvp_mean(ii,jj) = mean(vp(:,ii).*vp(:,jj),1);
    end
end

    end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAVE DATA

cd('C:\Users\rneuhausler\Everglades Research\filtered_Nov2013data')

if isProfiler
    t = t_clean;
    v = v_clean;
    pct_rm = pct_rm_clean;
    mint = min(t_clean{1});
    for c = 1:nCells, mint = min(mint,min(t_clean{c})); end
end

% print out percent removed
pct_rm

if ~exist('msg') % don't save any data if we produced an error message
if ~isempty(t)
    if ~isProfiler
        figure, plot(t, v),% xlim([0,min(10, max(t))])
        mint = t(1);
    end
    if calc_turbulence == 'y'
        save(to_save, 'v', 't', 'vpvp_mean','pct_rm')
    else
        save(to_save, 'v', 't','pct_rm')
    end
    fprintf('%s%d', 'First valid time(min)is ',mint)
else
    figure, plot(counter, v)
    if calc_turbulence == 'y'
        save(to_save, 'v', 'counter', 'vpvp_mean','pct_rm')
    else
        save(to_save, 'v', 'counter','pct_rm')
    end
    fprintf('%s%d', 'First valid count is ', counter(1))
end
end

end