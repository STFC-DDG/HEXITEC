%HEXITEC_Calibration_Generation_Dec_2020.m

%M. Veale 18/12/2020

%This script uses the calculated calibration values for a given module and
%then uses these the calibrate the response of each pixel before doing some
%basic analysis of the detector performance at a pixel level.

clear all
close all

%Open the files for analysis. NOTE: more than one file must be selected.

[fileName, pathName] = uigetfile('*.bin', 'Select the HEXITEC BIN files for analysis');
if isequal(fileName,0) || isequal(pathName,0)
    disp('ERROR - Check the number of files selected!')
    return
end

%Set the noise multiplier for threshold calculation

noise_multiplier = 6.5; %4 for CdTe @ 500V
Hit_Pixel_Threshold = 200;
highThresh = 8000;
binWidth = 10;
bins =0:binWidth:8000;

%% Noise Threshold Calculation

%The first step is to determine the optimal noise threshold per pixel. This
%is done by playing back the first chunk of the file and looking at the
%events that ocur at low signal values (<100 ADU). The peak of this
%distribution is then used to define the noise threshold.

tic

%This code locates the file and calculates it's size.

fileaddy = fullfile([pathName, char(fileName)]);
fid=fopen(fileaddy);
fseek(fid,0,'eof');
nFrames = floor(ftell(fid)/(80*80*2));
frewind(fid);
flag = 1;

%This preallocates the arrays for the threshold values.

noise_histo = zeros(80,80,100);
threshold = zeros(80,80);

%This set of loops completes a noise analysis for the first 50k frames
%of each of the selected files.

for fNum = 1:50000
    
    %This reads the current frame of data using the function written
    %by M. Wilson read80x80Frame.m.
    
    frame = read80x80Frame(fid,flag);
    
    %This loop inspects the value in each frame assuming it isn't a
    %real signal (<100ADU) and ats it to an array of noise histograms.
    
    for x = 1:80
        for y = 1:80
            if frame(x,y) > 0 && frame(x,y) < 100
                noise = frame(x,y);
                noise_histo(x,y,noise) = noise_histo(x,y,noise) + 1;
            end
        end
    end
end

%Having assembled the noise histograms for each pixel of this file an
%analysis of the ideal threshold position is made from a fit to the
%histogram data.

for x = 1:80
    for y = 1:80
        noise_plot = squeeze(noise_histo(x,y,:));
        noise_max = max(noise_plot);
        half_noise = round(noise_max/2);
        threshold(x,y) = find(noise_plot < half_noise,1,'first')*noise_multiplier;
    end
end

%This closes the current file which is open. This is important to avoid
%memory issues with Matlab.

fclose(fid);

toc

%% Data Analysis

tic

%Having calculated a low energy threshold for each pixel the analysis of
%the real data files can be completed. In this section of code charge
%sharing discrimination will be completed which requires the creation of
%some new frames and threshold matrices.

frame =zeros(80,80);
frame2 =zeros(80,80);
flag = 1;

%This block opens the file, calculates the number of frames and then
%rewinds to the beginning of the file ready for analysis.

fileaddy = fullfile([pathName, char(fileName)]);
fid=fopen(fileaddy);

%This preallocates the outputs of the script

RAW_Spec = zeros(80,80,length(bins));
CSD_Spec = zeros(80,80,length(bins));
RAW_Global_Spec = zeros(length(bins),1);
CSD_Global_Spec = zeros(length(bins),1);
multiplicity = zeros(80,80,6400);
No_Pixels = zeros(nFrames,1);
No_Clusters = zeros(nFrames,1);

%This loop now goes through each frame applying the energy calibration and
%thresholding based on the input files. Data is also processed for charge
%sharing.

for fNum = 1:nFrames
    
    frame = read80x80Frame(fid,flag);
    
    for x = 1:80
        for y = 1:80
            frame2(x,y) = frame(x,y);
            if frame(x,y) < threshold(x,y)
                frame2(x,y) = 0;
            end
            if frame2(x,y) <= 0
                frame2(x,y) = 0;
            else
                RAW_bin = round(frame2(x,y)/binWidth);
                if RAW_bin > 0 && RAW_bin < length(bins)
                    RAW_Global_Spec(RAW_bin,1) = RAW_Global_Spec(RAW_bin,1) + 1;
                    RAW_Spec(x,y,RAW_bin) = RAW_Spec(x,y,RAW_bin) + 1;
                end
                
            end
        end
    end
    
    No_Pixels(fNum,1) = length(nonzeros(frame2));
    
    frame3 = frame2;
    frame3(frame3>0) = 1;
    
    %This set of loops searches within the real data for hit pixels.
    
    Cluster_Search = bwconncomp(frame3);
    Clusters = labelmatrix(Cluster_Search);
    No_Clusters(fNum,1) = max(max(Clusters));
    
    for cl = 1:No_Clusters(fNum,1)
        [X_Hit,Y_Hit] = find(Clusters==cl);
        N_in_Cluster = length(X_Hit);
        
        Hits = zeros(length(X_Hit),1);
        
        for ncl = 1:length(X_Hit)
            Hits(ncl,1) = frame2(X_Hit(ncl,1),Y_Hit(ncl,1));
        end
        
        Max_Pixel = find(Hits==max(Hits));
        
        multiplicity(X_Hit(Max_Pixel,1),Y_Hit(Max_Pixel,1),N_in_Cluster) = multiplicity(X_Hit(Max_Pixel,1),Y_Hit(Max_Pixel,1),N_in_Cluster) + 1;
        
        if N_in_Cluster == 1
            CSD_bin = round(frame2(X_Hit,Y_Hit)/binWidth);
            if CSD_bin < length(bins) && CSD_bin > 0
                CSD_Spec(X_Hit,Y_Hit,CSD_bin) = CSD_Spec(X_Hit,Y_Hit,CSD_bin) + 1;
                CSD_Global_Spec(CSD_bin,1) = CSD_Global_Spec(CSD_bin,1) + 1;
            end
        end
         
    end
end

fclose(fid);
toc

%% Calibration Calculation

%In this section of the code the calculated histograms are used in order to
%complete a linear calibration per pixel and extract the fit coefficents.

%It is best to use the cleanest data set for the calibration of the
%detector, typically this is the discriminated data. Depending on the flux
%on the detector you may wish to use either the 5x5 or 3x3 data. Simply
%adjust the line of code below to switch between the two.

% Set up the matrices for peak amplitudes and centroids.

Peak_Amp = zeros(80,80,4);
Peak_Cent = zeros(80,80,4);
Rsq = zeros(80,80);
Max_Amp = zeros(80,80);
Max_Amp_Pos = zeros(80,80,2);

% These variables are used for the fitting. The first is the number of
% points used in the interpolation of individual spectra, i.e. the number
% of interpolation points between each data point. The other variable set
% the search window for the Am-241 primary photo peak.

Ninterp = 10;
Am_low = 1500; %in ADUs not bins!
Am_high = 2500; %in ADUs not bins!


% This loop goes through each spectrum, interpolates it and then searches
% for the Am-241 photo peak and identifies the peak centroid.

m_test = zeros(80,80);
c_test = zeros(80,80);

for x = 1:80
    for y = 1:80
        Test_Cal = zeros(2,2);
        spec_now = squeeze(CSD_Spec(x,y,:));
        spec_int = interp(spec_now,Ninterp);
        Peak_Amp(x,y,1) = max(spec_int(Am_low:Am_high,1));
        Max_Amp_Pos(x,y,1) = (Am_low+find(spec_int(Am_low:Am_high,1)==Peak_Amp(x,y,1),1,'first'));
        Max_Amp(x,y) = max(spec_int(200:1000));
        Max_Amp_Pos(x,y,2) = (200+find(spec_int(200:1000,1)==Max_Amp(x,y),1,'first'));
        Test_Cal(:,1) = [59.54 17.75];
        Test_Cal(1,2) = Max_Amp_Pos(x,y,1);
        Test_Cal(2,2) = Max_Amp_Pos(x,y,2);
        [coeffs] = polyfit(Test_Cal(:,2),Test_Cal(:,1),1);
        m_test(x,y) = coeffs(1,1);
        c_test(x,y) = coeffs(1,2);
    end
end

% Having calculated the Am-241 peak position it's possible to make a first
% rough estimate of the calibration of each pixel. This is then used to
% guess the position of the three other peaks and then searches for the
% actual position in a window around the guess.

% Am_Cent = squeeze(Peak_Cent(:,:,1));
% m_rough = E_Peak(1,1)./Am_Cent;
% Peak_Window_Low(:,1) = [100 25 25 25]; %in ADUs not bins!
% Peak_Window_High(:,1) = [100 25 25 25]; %in ADUs not bins!

E_Peak(:,1) = [59.54 13.94 17.75 26.34];
E_window = 1; %keV
Peak_Window_Low = zeros(4,1);
Peak_Window_High = zeros(4,1);

for E = 1:4
    Peak_Window_Low(E,1) = round((E_Peak(E,1)-E_window-c_test(x,y))/m_test(x,y));
    Peak_Window_High(E,1) = round((E_Peak(E,1)+E_window-c_test(x,y))/m_test(x,y));
end

m = zeros(80,80);
c = zeros(80,80);

for x= 1:80
    for y = 1:80
        
        spec_now = squeeze(CSD_Spec(x,y,:));
        spec_int = interp(spec_now,Ninterp);
        
        if min(Peak_Window_Low) > -1
            
            for z = 1:4
                Peak_min = Peak_Window_Low(z,1);
                Peak_max = Peak_Window_High(z,1);
                Peak_Amp(x,y,z) = max(spec_int(Peak_min:Peak_max,1));
                Peak_Cent(x,y,z) = (Peak_min+find(spec_int(Peak_min:Peak_max,1)==Peak_Amp(x,y,z),1,'first'));
            end
            
            %Extracts the calculations for the peak centroids in order to
            %perform a linear fit to the data.
            
            ADU_Peak(:,1) = squeeze(Peak_Cent(x,y,:));
            [coeffs] = polyfit(ADU_Peak,E_Peak,1);
            R = corrcoef(ADU_Peak,E_Peak);
            Rsq(x,y) = R(1,2).^2;
            
%             if Rsq(x,y) < 0.98
%                 
%                 figure (1)
%                 scatter(ADU_Peak,E_Peak)
%                 figure (2)
%                 plot(spec_int(1:2500))
%                 display(ADU_Peak)
%                 pause
%                 
%             end
            
            m(x,y) = coeffs(1,1);
            c(x,y) = coeffs(1,2);
            
            %         if x == 32 && y == 62
            %             pause
            %         end
        end
    end
end

%% Data Visulisation

figure (1)
imagesc(m*1000)
colorMap = jet(256);
colormap(colorMap)
set(gca,'FontSize',18,'FontWeight','bold','fontname','Arial')
caxis([25 30])
d = colorbar;
d.Label.String = 'Gradient (eV ADU^-^1)';
d.Label.FontSize = 22;
xlabel('X (Pix)','FontSize', 22,'FontWeight','bold')
ylabel('Y (Pix)','FontSize', 22,'FontWeight','bold')
axis square

figure (2)
imagesc(c)
colorMap = jet(256);
colormap(colorMap)
set(gca,'FontSize',18,'FontWeight','bold','fontname','Arial')
caxis([0 4])
d = colorbar;
d.Label.String = 'Intercept (keV)';
d.Label.FontSize = 22;
xlabel('X (Pix)','FontSize', 22,'FontWeight','bold')
ylabel('Y (Pix)','FontSize', 22,'FontWeight','bold')
axis square

E_thresh = zeros(80,80);

for x = 1:80
    for y = 1:80
        val = threshold(x,y,1);
        E_val = (m(x,y)*val)+c(x,y);
        E_thresh(x,y) = E_val;
    end
end

figure (3)
imagesc(threshold)
colorMap = jet(256);
colormap(colorMap)
set(gca,'FontSize',18,'FontWeight','bold','fontname','Arial')
caxis([0 80])
d = colorbar;
d.Label.String = 'Threshold (ADU)';
d.Label.FontSize = 22;
xlabel('X (Pix)','FontSize', 22,'FontWeight','bold')
ylabel('Y (Pix)','FontSize', 22,'FontWeight','bold')
axis square

figure (4)
imagesc(Rsq)
colorMap = jet(256);
colormap(colorMap)
set(gca,'FontSize',18,'FontWeight','bold','fontname','Arial')
caxis([0.9998 1.0])
d = colorbar;
d.Label.String = 'R^2';
d.Label.FontSize = 22;
xlabel('X (Pix)','FontSize', 22,'FontWeight','bold')
ylabel('Y (Pix)','FontSize', 22,'FontWeight','bold')
axis square

figure (5)
plot(No_Pixels,'Color', 'k')
set(gca,'FontSize',18,'FontWeight','bold','fontname','Arial')
xlabel('Frame No.','FontSize', 22,'FontWeight','bold')
ylabel('Hits','FontSize', 22,'FontWeight','bold')
