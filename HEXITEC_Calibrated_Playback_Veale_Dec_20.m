%Dec20_HEXITEC_Calibrated_Playback.m

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

%Load in the linear calibration coefficent arrays

load('Veale_m.mat');
load('Veale_c.mat');
load('Veale_threshold.mat');

%The high energy threshold sets the maximum bin size. Typically a bin width
%of 0.25 keV is used for CdTe and CdZnTe detectors, 0.5 keV is more
%suitable for GaAs detectors.

highThresh = 250;
binWidth = 0.25;
bins =0:binWidth:highThresh;
Hit_Pixel_Threshold = 200;

%% Data Analysis

%This prealocates the frames that are used in the analysis.

frame =zeros(80,80);
frame2 =zeros(80,80);
frame3 =zeros(80,80);
flag = 1;

tic

%This block opens the file, calculates the number of frames and then
%rewinds to the beginning of the file ready for analysis.

fileaddy = fullfile([pathName, char(fileName)]);
fid=fopen(fileaddy);
fseek(fid,0,'eof');
nFrames = floor(ftell(fid)/(80*80*2));
frewind(fid);

%The analysis output arrays are preallocated.

RAW_Spec = zeros(80,80,length(bins));
CSD_Spec = zeros(80,80,length(bins));
CSA_Spec = zeros(80,80,length(bins)*2,3);

RAW_Global_Spec = zeros(length(bins),1);
CSD_Global_Spec = zeros(length(bins),1);
CSA_Global_Spec = zeros(length(bins)*2,3);

multiplicity = zeros(80,80,6400);
No_Pixels = zeros(nFrames,1);
CS_Dist = zeros(length(bins)*2,2000);
CS_Dist_Multi = zeros(length(bins)*4,length(bins)*4,3);

%This loop now goes through each frame applying the energy calibration and
%thresholding based on the input files. Data is also processed for charge
%sharing.

for fNum = 1:nFrames
    
    frame = read80x80Frame(fid,flag);
    
    for x = 1:80
        for y = 1:80
            frame2(x,y) = (frame(x,y)*m(x,y))+c(x,y);
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
    No_Clusters = max(max(Clusters));
    
    for cl = 1:No_Clusters
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
        elseif N_in_Cluster == 2
            E1 = frame2(X_Hit(1,1),Y_Hit(1,1));
            E1_bin = round(E1/binWidth);
            
            E2 = frame2(X_Hit(2,1),Y_Hit(2,1));
            E2_bin = round(E2/binWidth);
            
            ET = E1 + E2;
            CSA_bin = round(ET/binWidth);
            
            Ratio = round((((E1-E2)/ET)+1)/0.001);
            
            if CSA_bin < length(bins)*2 && CSA_bin > 0
                if E1 > E2
                    CSA_Spec(X_Hit(1,1),Y_Hit(1,1),CSA_bin,1) = CSA_Spec(X_Hit(1,1),Y_Hit(1,1),CSA_bin,1) + 1;
                    CS_Dist_Multi(E1_bin,E2_bin,1) = CS_Dist_Multi(E1_bin,E2_bin,1) + 1;
                else
                    CSA_Spec(X_Hit(2,1),Y_Hit(2,1),CSA_bin,1) = CSA_Spec(X_Hit(2,1),Y_Hit(2,1),CSA_bin,1) + 1;
                    CS_Dist_Multi(E2_bin,E1_bin,1) = CS_Dist_Multi(E2_bin,E1_bin,1) + 1;
                end
                CSA_Global_Spec(CSA_bin,1) = CSA_Global_Spec(CSA_bin,1) + 1;
                CS_Dist(CSA_bin,Ratio) = CS_Dist(CSA_bin,Ratio) + 1;
            end
        elseif N_in_Cluster == 3
            E_Max = Hits(Max_Pixel,1);
            E_Neighbours = sum(Hits)-E_Max;
            ET = sum(Hits);
            Em = round(E_Max/binWidth);
            En = round(E_Neighbours/binWidth);
            CSA_bin = round(ET/binWidth);
            if CSA_bin < length(bins)*2 && CSA_bin > 0
                CS_Dist_Multi(Em,En,2) = CS_Dist_Multi(Em,En,2) + 1;
                CSA_Global_Spec(CSA_bin,2) = CSA_Global_Spec(CSA_bin,2) + 1;
                CSA_Spec(X_Hit(Max_Pixel,1),Y_Hit(Max_Pixel,1),CSA_bin,2) = CSA_Spec(X_Hit(Max_Pixel,1),Y_Hit(Max_Pixel,1),CSA_bin,2) + 1;
            end
        elseif N_in_Cluster == 4
            E_Max = Hits(Max_Pixel,1);
            E_Neighbours = sum(Hits)-E_Max;
            ET = sum(Hits);
            Em = round(E_Max/binWidth);
            En = round(E_Neighbours/binWidth);
            CSA_bin = round(ET/binWidth);
            if CSA_bin < length(bins)*2 && CSA_bin > 0
                CS_Dist_Multi(Em,En,3) = CS_Dist_Multi(Em,En,3) + 1;
                CSA_Global_Spec(CSA_bin,3) = CSA_Global_Spec(CSA_bin,3) + 1;
                CSA_Spec(X_Hit(Max_Pixel,1),Y_Hit(Max_Pixel,1),CSA_bin,3) = CSA_Spec(X_Hit(Max_Pixel,1),Y_Hit(Max_Pixel,1),CSA_bin,3) + 1;
            end
        end
         
    end
end

fclose(fid);
toc

%% Multiplicity Analysis

%This section simply takes the multiplicity data per pixel and calculates
%an average value per pixel.

pixel_sharing = zeros(80,80);
Multiplicity_Histo = zeros(6400,1);

for x = 1:80
    for y = 1:80
        
        multi_now(:,1) = squeeze(multiplicity(x,y,:));
        multi_counts = sum(multi_now);
        total_multi = 0;
        pixel_sharing(x,y) = (1-(multi_now(1,1)/multi_counts))*100;
        
        for z = 1:length(multi_now)
            total_multi = total_multi + (z*multi_now(z,1));
        end
        
        Multiplicity_Histo = Multiplicity_Histo + multi_now;
    end
end

figure (1)
imagesc(pixel_sharing)
colorMap = jet(256);
colormap(colorMap)
set(gca,'FontSize',18,'FontWeight','bold','fontname','times')
caxis([10 70])
d = colorbar;
% d.Label.String = 'Charge Sharing Percentage (%)';
% d.Label.FontSize = 22;
title('Charge Sharing (%)','FontSize', 22,'FontWeight','bold')
xlabel('X (Pix)','FontSize', 22,'FontWeight','bold')
ylabel('Y (Pix)','FontSize', 22,'FontWeight','bold')
axis square

%% SPECTROSCOPIC PERFORMANCE

%In this section of the code the energy resolution of the detector and the
%pixel yields are calculated.

Ninterp = 10;
E_Peak = 59.54;   %keV
E_Window = 2;   %keV

%Sets up the peak matrices

Peak_Counts = zeros(80,80);
Peak_Amp = zeros(80,80);
Peak_Cent = zeros(80,80);
Peak_LWHM = zeros(80,80);
Peak_UWHM = zeros(80,80);
Peak_FWHM = zeros(80,80);

for x = 1:80
    for y = 1:80       
        E_low = round(((E_Peak-E_Window)/binWidth))*Ninterp;
        E_high = round(((E_Peak+E_Window+10)/binWidth))*Ninterp;

        spec = squeeze(CSD_Spec(x,y,:));
        spec_int = interp(spec,Ninterp);

        Peak_Counts(x,y) = sum(spec_int(E_low:E_high,1));
        
        if Peak_Counts(x,y) > 100
            Peak_Amp(x,y) = max(spec_int(E_low:E_high));
            Peak_Cent(x,y) = (E_low+find(spec_int(E_low:E_high,1)==Peak_Amp(x,y),1,'first'));
            Peak_LWHM(x,y) = (E_low+find(spec_int(E_low:Peak_Cent(x,y),1)>Peak_Amp(x,y)/2,1,'first'));
            Peak_UWHM(x,y) = (Peak_Cent(x,y)+find(spec_int(Peak_Cent(x,y):E_high,1)<Peak_Amp(x,y)/2,1,'first'));
            FWHM = (Peak_UWHM(x,y) - Peak_LWHM(x,y));
            if FWHM > 0
                Peak_FWHM(x,y) = ((Peak_UWHM(x,y) - Peak_LWHM(x,y))/Ninterp)*binWidth;
            end
        end
        
    end
end

%% Data Visualisation

figure (2)
imagesc(Peak_FWHM)
colorMap = jet(256);
colormap(colorMap)
set(gca,'FontSize',18,'FontWeight','bold','fontname','times')
caxis([0 2])
d = colorbar;
d.Label.String = '59.54keV FWHM (keV)';
d.Label.FontSize = 22;
%title('FWHM (keV)','FontSize', 22,'FontWeight','bold')
xlabel('X (Pix)','FontSize', 22,'FontWeight','bold')
ylabel('Y (Pix)','FontSize', 22,'FontWeight','bold')
axis square

figure (3)
imagesc(Peak_Counts)
colorMap = jet(256);
colormap(colorMap)
set(gca,'FontSize',18,'FontWeight','bold','fontname','times')
caxis([0 30000])
d = colorbar;
d.Label.String = '59.54keV Counts';
d.Label.FontSize = 22;
xlabel('X (Pix)','FontSize', 22,'FontWeight','bold')
ylabel('Y (Pix)','FontSize', 22,'FontWeight','bold')
axis square

figure (4)
imagesc(log10(CS_Dist(1:500,:)))
colorMap = jet(256);
colormap(colorMap)
set(gca,'FontSize',18,'FontWeight','bold','fontname','times','YDir','normal')
caxis([0 3])
d = colorbar;
d.Label.String = 'Counts';
d.Label.FontSize = 22;
xlabel('(E_{1}-E_{2})/E_{T}','FontSize', 22,'FontWeight','bold')
ylabel('E_{T} (keV)','FontSize', 22,'FontWeight','bold')
xticks([1 500 1000 1500 2000])
xticklabels([-1.0 -0.5 0 0.5 1.0])
yticks([1 100 200 300 400 500])
yticklabels([0 25 50 75 100 125])
axis square

figure (5)
imagesc(log(CS_Dist_Multi(1:300,1:300,1)))
colorMap = jet(256);
colormap(colorMap)
set(gca,'FontSize',18,'FontWeight','bold','fontname','times','YDir','normal')
caxis([1 10])
d = colorbar;
d.Label.String = 'Log_{10}(Counts)';
d.Label.FontSize = 22;
xlabel('(E_{1}-E_{2})/E_{T}','FontSize', 22,'FontWeight','bold')
ylabel('E_{T} (keV)','FontSize', 22,'FontWeight','bold')
axis square