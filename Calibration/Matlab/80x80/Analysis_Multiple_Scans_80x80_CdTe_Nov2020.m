%Script to process the ISIS Gd Data
% Take raw - reeduced format data
% sort to spectra per pixel
% save as hxt file

clear all
close all

%% Accessing files to analyse
fileString = 'J:\NewCdTe'
cd (fileString)
filesList = dir('1606*'); %Detector stub val goes here
% filesList = filesList(5:8);

%Change depending on detector properties
noiseThreshold = 80;

%%
%loops through the directories containing each bin file for detector
for sensorNum=1:length(filesList)
    cd (fileString)
    close all
    clearvars -except filesList sensorNum noiseThreshold fileString
    addpath 'H:\Documents\Hexitec\Matlab'
    
    %opening of directory containing bin file
    sensorName = filesList(sensorNum).name;
    pathName = filesList(sensorNum).folder;

    sensorPath = strcat(pathName,'\',sensorName);
    
%     sensorPath = strcat(sensorPath,'\','500V');
    
    cd (sensorPath)
    
    %DEFAULTS TO SEARCHING FOR 750v FILE, IF UNAVAILABLE SEARCHES FOR 600V,
    %IF THAT'S NOT AVAILIBLE SEARCHES FOR 500V
    file750V = dir('*400V*.bin');
    
    if isempty(file750V)==1
        file600V = dir('*600V*.bin');
        file750V = file600V;
        if isempty(file600V)==1
            file500V = dir('*500V*.bin');
            file750V = file500V;
        end
    end

    fileSizes = zeros(length(file750V),1);
    if length(file750V)>1
        for z=1:length(file750V)
            s=dir(strcat(file750V(z).folder,'\',file750V(z).name));
            fileSizes(z) = s.bytes;
        end
    end

    bigFile = max(fileSizes);
    fullFile = find(fileSizes == bigFile, 1, 'first');

    currentFile = file750V(fullFile);
    filePath = strcat(currentFile.folder,'\',currentFile.name);
    
    %this is the file that now goes on to be analysed and frames processed
    %from.
    currentFile
    %%
    fid=fopen(filePath);

    fseek(fid,0,'eof');
    nFrames = floor(ftell(fid)/(80*80));
    frewind(fid);
    flag = 1;

    %preallocation of variables
    frame =zeros(80,80);
    
    binMin=1;
    binMax=4000;
    binWidth=5;
    edges=binMin:binWidth:binMax;
    globalHist=zeros(length(edges),1);
    globalHistCSD=zeros(length(edges),1);
    globalHistCSA=zeros(length(edges),1);
    pixelHist=zeros(length(edges),80,80);
    pixelHistCSA=zeros(length(edges),80,80);
    pixelHistCSD=zeros(length(edges),80,80);
    noiseHist=zeros(length(edges),80,80);
    noiseMultiplier=8;
    noiseMatrix=zeros(80,80);

    %first quick loop of the files..... I'm not entirely sure why? For
    %finding the noise?
    for fNum = 1:10000%(nFrames/2)-1
        frame = read80x80Frame(fid,flag);

        %bins the frame values
        binnedFrame = discretize(frame,edges);

        %puts bins in individual pixel histograms
        [x,y]=find(binnedFrame>0);
        for z=1:length(x)
            pixelHist(binnedFrame(x(z),y(z)))=pixelHist(binnedFrame(x(z),y(z)),x(z),y(z))+1;
            if binnedFrame(x(z),y(z))>0 && binnedFrame(x(z),y(z))<100
                noiseHist(binnedFrame(x(z),y(z)),x(z),y(z))=noiseHist(binnedFrame(x(z),y(z)),x(z),y(z))+1;
            end
        end

        %puts bins in global histogram
        binnedFrame=binnedFrame(binnedFrame>0);
        for z=1:length(binnedFrame)
            globalHist(binnedFrame(z))=globalHist(binnedFrame(z))+1;
        end
    end

    %finds the location of the half noise and defines a noise threshold based
    %on pre defined noiseMultiplier value
    noiseHalf=(squeeze(max(noiseHist,[],1)))/2;
    [x,y] = find(noiseHalf>0);
    for z=1:length(x)
        noiseMatrix(x(z),y(z))=find(noiseHist(:,x(z),y(z))<round(noiseHalf(x(z),y(z))),1,'first')*noiseMultiplier;
    end

    fclose(fid);

    if sum(any(noiseMatrix==0))>0
        [x,y]=find(noiseMatrix==0);
        meanNoise = round(mean(noiseMatrix(noiseMatrix>0)));
        for z=1:length(x)
            noiseMatrix(x(z),y(z))=meanNoise;
        end
    end

    %% This is the analysis section following thresholding with CSD and CSA - Rhian's Method!!
    tic
    fid=fopen(filePath);

    fseek(fid,0,'eof');
    nFrames = floor(ftell(fid)/(80*80));
    frewind(fid);
    flag = 1;
    
    globalHist=zeros(length(edges),1);
    pixelHist=zeros(length(edges),80,80);
    
    %now properly reads all the frames of data!
    for fNum = 1:(nFrames/4)-1
        frame = read80x80Frame(fid,flag);
        frameCSD = zeros(80,80);
        frameCSA = zeros(80,80);
        
        %bins the frame values
        binnedFrame = discretize(frame,edges);

        %puts bins in individual pixel histograms
        [x,y]=find(binnedFrame>0);
        for z=1:length(x)
            pixelHist(binnedFrame(x(z),y(z)),x(z),y(z))=pixelHist(binnedFrame(x(z),y(z)),x(z),y(z))+1;
            if binnedFrame(x(z),y(z))>0 && binnedFrame(x(z),y(z))<100
                noiseHist(binnedFrame(x(z),y(z)),x(z),y(z))=noiseHist(binnedFrame(x(z),y(z)),x(z),y(z))+1;
            end
        end

        %puts bins in global histogram
        binnedFrame=binnedFrame(binnedFrame>0);
        for z=1:length(binnedFrame)
            globalHist(binnedFrame(z))=globalHist(binnedFrame(z))+1;
        end

        %CSD and CSA
        %double pass of clusters identification to stop over counting in
        %CSA - don't believe this is actually necessary now...but if it
        %ain't broke don't fix it!!
        clusterMaskAll = frame>noiseThreshold;%noiseMatrix;
        clusterMaskLow = immultiply((frame>noiseMatrix),(frame<1000));
        clusterMaskHigh = (frame>1000);

        clusterMask=clusterMaskAll;

        frameThresh = immultiply(clusterMask,frame);
        frameCSD = frameThresh;
        frameCSA = frameThresh;
        
        %this function labels the clusters numericall, we loop through each
        %cluster in the next group.
        clustersImage = bwlabel(clusterMask);
        numClusters = max(max(clustersImage));
        
        %big ass loop below essentally goes through all the clusters
        %identified in turn. CSD chucks them if they're > 1 pixel in size.
        %CSA chucks if they're greater than 4 in size. CSA also chucks if
        %difference between 2 points in cluster >500 ADU (because probably
        %not charge sharing at that level. 
        for currentCluster=1:numClusters
            [x,y]=find(clustersImage==currentCluster);
            if length(x)>1
                for z=1:length(x)
                    frameCSD(x(z),y(z))=0;
                end
                totalCharge=sum(unique(frameThresh(x,y)));
                largeFrac=max(unique(frameThresh(x,y)));
                smallFrac=min(unique(frameThresh(x,y)));
                [i,j]=find(frameThresh==largeFrac);
                if length(i)>1
                    for z=1:length(i)
                    testElements = find(x==i(z));
                    if isempty(testElements)==0
                        coNum = z;
                        break
                    end
                    end                
                for z=1:length(x)
                    frameCSA(x(z),y(z))=0;
                end   
                frameCSA(i(coNum),j(coNum))=totalCharge;
                else
                    for z=1:length(x)
                        frameCSA(x(z),y(z))=0;
                    end
                    frameCSA(i,j)=totalCharge;
                end
                if (largeFrac-smallFrac)>500
                    for z=1:length(x)
                        frameCSA(x(z),y(z))=0;
                    end
                end
            end
            if length(x)>4
                for z=1:length(x)
                    frameCSA(x(z),y(z))=0;
                end
            end
        end

        %bins the frame values for CSD
        binnedFrameCSD = discretize(frameCSD,edges);

        %puts CSD bins in individual pixel histograms
        [x,y]=find(binnedFrameCSD>0);
        for z=1:length(x)
            pixelHistCSD(binnedFrameCSD(x(z),y(z)),x(z),y(z))=pixelHistCSD(binnedFrameCSD(x(z),y(z)),x(z),y(z))+1;
        end

        %puts CSD bins in global histogram
        binnedFrameCSD=binnedFrameCSD(binnedFrameCSD>0);
        for z=1:length(binnedFrameCSD)
            globalHistCSD(binnedFrameCSD(z))=globalHistCSD(binnedFrameCSD(z))+1;
        end

        %bins the frame values for CSA
        binnedFrameCSA = discretize(frameCSA,edges);

        %puts CSD bins in individual pixel histograms
        [x,y]=find(binnedFrameCSA>0);
        for z=1:length(x)
            pixelHistCSA(binnedFrameCSA(x(z),y(z)),x(z),y(z))=pixelHistCSA(binnedFrameCSA(x(z),y(z)),x(z),y(z))+1;
        end

        %puts CSD bins in global histogram
        binnedFrameCSA=binnedFrameCSA(binnedFrameCSA>0);
        for z=1:length(binnedFrameCSA)
            globalHistCSA(binnedFrameCSA(z))=globalHistCSA(binnedFrameCSA(z))+1;
        end

    end
    fclose(fid)
    toc

    [pks,locs,widths,prom] = findpeaks(globalHistCSD);

    limit=800;

    figure
    globalEdges = edges;
    plot(globalEdges,globalHistCSD)

disp('if this is first run through EXIT CODE NOW and use command: plot(globalEdges,globalHistCSD) to identify upper and lower energy peak limits for calibration')
% pause

    %% Calibration
    %Calibration on CSD

    FitSpec = pixelHistCSD;
    peakAmplitudes = zeros(80,80,3);
    peakCentroids = zeros(80,80,3);
    m = zeros(80,80);
    c = zeros(80,80);
    image60keV = zeros(80,80);

    AmLow = 1500; %ADU!
    AmHigh = 2500; %ADU!

    %interpolation of data required later to allow use of values in ADU and
    %prevent overlapping peaks. nInterp == binWidth to convert back to ADU.
    nInterp = binWidth;

    energyPeaks = [59.54 13.94 17.75];
    searchLow = [1500 350 500];
    searchHigh = [2000 499 650];

    y=1;
    x=0;
    for z = 1:(80*80)
        x=x+1;

        currentPixel = interp(squeeze(FitSpec(:,x,y)),nInterp);
        peakAmplitudes(x,y,1) = max(currentPixel(AmLow:AmHigh));
        peakCentroids(x,y,1) = AmLow + (find(currentPixel(AmLow:AmHigh)==peakAmplitudes(x,y,1),1,'first'));
        halfWidth = 100;
        image60keV(x,y) = sum(currentPixel(peakCentroids(x,y,1)-halfWidth:peakCentroids(x,y,1)+halfWidth));

        if x==80
            x=0;
            y=y+1;
        end
    end

    mRough = squeeze(peakCentroids(:,:,1))/energyPeaks(1);
    lowWindow = 50;
    highWindow = 50;

    y=1;
    x=0;
    for z = 1:(80*80)
        x=x+1;

        currentPixel = interp(squeeze(FitSpec(:,x,y)),nInterp);
        for peakNum=2:3        
            peakAmplitudes(x,y,peakNum)=max(currentPixel(searchLow(peakNum):searchHigh(peakNum)));
            peakCentroids(x,y,peakNum)=searchLow(peakNum)+ (find(currentPixel(searchLow(peakNum):searchHigh(peakNum))==peakAmplitudes(x,y,peakNum),1,'first'));
        end

        fitData = squeeze(peakCentroids(x,y,:));
        fitData = fitData';
        [coeffs] = polyfit(fitData,energyPeaks,1);
        m(x,y) = coeffs(1,1);
        c(x,y) = coeffs(1,2);

        if x==80
            x=0;
            y=y+1;
        end
    end

    %% Data Viulisation
    figureFile = strcat('Figures_400V_', sensorName);
    
    %makes a directory to save all the figures in
    mkdir(figureFile)
    cd(figureFile)

    fig = figure (1) ;
    imagesc(m) 
    colorMap = jet(256); 
    colormap(colorMap) 
    set(gca,'FontSize',18,'fontname','times') 
    caxis([min(min(m)) max(max(m))]) 
    d = colorbar; 
    d.Label.String = 'Gradient (keV ADU-1)'; 
    d.Label.FontSize = 18; 
    xlabel('X (Pix)','FontSize', 18) 
    ylabel('Y (Pix)','FontSize', 18) 
    axis square 
    saveas(fig,'Gradients.fig')

    fig = figure (2) ;
    imagesc(c) 
    colorMap = jet(256); 
    colormap(colorMap) 
    set(gca,'FontSize',18,'fontname','times') 
    caxis([min(min(c)) max(max(c))]) 
    d = colorbar; 
    d.Label.String = 'Intercept (keV)'; 
    d.Label.FontSize = 18; 
    xlabel('X (Pix)','FontSize', 18) 
    ylabel('Y (Pix)','FontSize', 18) 
    axis square 
    saveas(fig,'Intercepts.fig')

    fig = figure (3) ;
    imagesc(noiseMatrix) 
    colorMap = jet(256); 
    colormap(colorMap) 
    set(gca,'FontSize',18,'fontname','times') 
    caxis([min(min(noiseMatrix)) max(max(noiseMatrix))]) 
    d = colorbar; 
    d.Label.String = 'Threshold (keV)'; 
    d.Label.FontSize = 18; 
    xlabel('X (Pix)','FontSize', 18) 
    ylabel('Y (Pix)','FontSize', 18) 
    axis square 
    saveas(fig,'Noise.fig')

    fig = figure (4);
    imagesc(image60keV) 
    colorMap = jet(256); 
    colormap(colorMap) 
    set(gca,'FontSize',18,'fontname','times') 
    caxis([0 10000]) 
    d = colorbar; 
    d.Label.String = '60 keV Counts'; 
    d.Label.FontSize = 18; 
    xlabel('X (Pix)','FontSize', 18) 
    ylabel('Y (Pix)','FontSize', 18) 
    axis square 
    saveas(fig,'Image60keV.fig')
    %% Histogrammed Data

    [n,edges]=histcounts(m,50);
    cumHist=cumsum(n);
    cumHist = (cumHist/max(cumHist))*100;
    fig = figure(5);
    yyaxis left
    a=histogram(m,edges, 'FaceColor', 'b', 'FaceAlpha', 0.9, 'EdgeColor', 'w', 'LineWidth', 0.8)
    set(gca,'FontSize',16,'fontname','times') 
    xlabel('Gradient (keV ADU^{-1})','FontSize', 16) 
    ylabel('# Pixels','FontSize', 16) 
    yyaxis right
    b=plot(edges(1:length(edges)-1),cumHist)
    ylabel('Cumulative Pixels as Percent','FontSize', 16) 
    b.LineWidth=3
    saveas(fig,'GradientHistogram.fig')

    [n,edges]=histcounts(c,50);
    cumHist=cumsum(n);
    cumHist = (cumHist/max(cumHist))*100;
    fig = figure(6);
    yyaxis left
    a=histogram(c,edges, 'FaceColor', 'b', 'FaceAlpha', 0.9, 'EdgeColor', 'w', 'LineWidth', 0.8)
    set(gca,'FontSize',16,'fontname','times') 
    xlabel('Intercept (keV)','FontSize', 16) 
    ylabel('# Pixels','FontSize', 16) 
    yyaxis right
    b=plot(edges(1:length(edges)-1),cumHist)
    ylabel('Cumulative Pixels as Percent','FontSize', 16) 
    b.LineWidth=3
    saveas(fig,'InterceptsHistogram.fig')

    [n,edges]=histcounts(noiseMatrix,50);
    cumHist=cumsum(n);
    cumHist = (cumHist/max(cumHist))*100;
    fig = figure(7);
    yyaxis left
    a=histogram(noiseMatrix,edges, 'FaceColor', 'b', 'FaceAlpha', 0.9, 'EdgeColor', 'w', 'LineWidth', 0.8)
    set(gca,'FontSize',16,'fontname','times') 
    xlabel('Threshold (keV)','FontSize', 16) 
    ylabel('# Pixels','FontSize', 16) 
    yyaxis right
    b=plot(edges(1:length(edges)-1),cumHist)
    ylabel('Cumulative Pixels as Percent','FontSize', 16) 
    b.LineWidth=3
    saveas(fig, 'ThresholdsHistogram.fig')

    %% Calibrated plots

    %individual pixel calibrations
    edges=binMin:binWidth:binMax;
    calibratedEdges = zeros(length(edges),80,80);
    x=0; y=1;
    for z = 1:(80*80)
        x=x+1;
        edgesCopy=edges;
        edgesCopy=(edgesCopy*m(x,y)) + c(x,y);
        calibratedEdges(:,x,y) = edgesCopy;    
        if x==80
            x=0;
            y=y+1;
        end
    end
    fig = figure (8)
    plot(squeeze(calibratedEdges(:,40,40)),squeeze(pixelHistCSD(:,40,40)), 'b')
    set(gca,'FontSize',18,'fontname','times') 
    xlabel('Energy (keV)','FontSize', 18)
    ylabel('Counts)','FontSize', 18)
    title('Single Pixel Spectra - CSD')
    box off
    saveas(fig,'CalibratedCSDplotPixel4040.fig')

    %global calibration
    globalEdges = squeeze(mean(mean(calibratedEdges,2),3));
    fig = figure (9)
    plot(globalEdges,globalHistCSD, 'b')
    set(gca,'FontSize',18,'fontname','times') 
    xlabel('Energy (keV)','FontSize', 18)
    ylabel('Counts)','FontSize', 18)
    title('Global Spectra - CSD')
    box off
    saveas(fig,'CalibratedGlobalCSDplot.fig')
    %% Peaks analysis
    peakCounts = zeros(80,80);
    peakAmp = zeros(80,80);
    peakCentADU = zeros(80,80);
    peakCentkeV = zeros(80,80);

    Epeak = 60; %keV
    searchWindow=5; %keV
    ELowOriginal = Epeak-searchWindow;
    EHighOriginal = Epeak+searchWindow;

    x=0; y=1;
    for z = 1:(80*80)
        x=x+1;
        
        ELow = find(globalEdges>ELowOriginal,1,'first');
        EHigh = find(globalEdges>EHighOriginal,1,'first');

        spec=squeeze(pixelHistCSD(ELow:EHigh,x,y));

        if isempty(spec)==1
            continue
        end
        specInterp = interp(spec,nInterp);

        peakCounts(x,y) = sum(specInterp);
        peakAmp(x,y)=max(specInterp);
        peakCentADU(x,y) = ((find(specInterp==(max(specInterp)),1,'first'))/nInterp)+ELow;
        peakCentkeV(x,y) = (peakCentADU(x,y)*m(x,y)) + c(x,y);


        if x==80
            x=0;
            y=y+1;
        end
    end

    fig = figure (10);
    imagesc(peakCentkeV)
    saveas(fig,'PeakCentroidsImage.fig');
    
    mainPeakCents = squeeze(peakCentroids(:,:,1));
    fig = figure(11);
    imagesc(mainPeakCents)
    colorMap = jet(256); 
    colormap(colorMap) 
    set(gca,'FontSize',18,'fontname','times') 
    caxis([min(min(mainPeakCents)) max(max(mainPeakCents))]) 
    d = colorbar; 
    d.Label.String = '60keV Peak Centroids (ADU)'; 
    d.Label.FontSize = 18; 
    xlabel('X (Pix)','FontSize', 18) 
    ylabel('Y (Pix)','FontSize', 18) 
    axis square 
    saveas(fig,'peakCents_60keV')
    
    %saving of useful matrices
    save('savepixelHistCSD.mat','pixelHistCSD')
    save('savepixelHistCSA.mat','pixelHistCSA')
    save('savepixelHist.mat','pixelHist')
    save('saveglobalHistCSD.mat','globalHistCSD')
    save('saveglobalHistCSA.mat','globalHistCSA')
    save('saveglobalHist.mat','globalHist')
    save('saveEdges.mat','edges')
    save('saveCalibratedEdges','calibratedEdges')
   
    %% Corrected FWHM Calculations
    x=0;
    peakHalf = zeros(80,80);
    FWHM = zeros(80,80);
    FWHMKeV = zeros(80,80);
    peakAmplitudes = zeros(80,80);
    peakCentroids = zeros(80,80);
    StepsX = 80;
    y=1;
    for z = 1:(StepsX*StepsX)
        x=x+1;

        currentSubPix = squeeze(pixelHistCSD(250:800,x,y));

        peakAmplitudes(x,y)=max(currentSubPix);
        peakCentroids(x,y)=find(currentSubPix==peakAmplitudes(x,y),1,'first');    
        peakHalf(x,y) =  peakAmplitudes(x,y)/2;
        FWHMHigh = peakCentroids(x,y) + (find(currentSubPix(peakCentroids(x,y):length(currentSubPix))< peakHalf(x,y),1,'first'));
        FWHMLow = find(currentSubPix>peakHalf(x,y), 1, 'first');
        if isempty(FWHMHigh)==1
            continue
        end
        FWHMKeV(x,y) = calibratedEdges(250+FWHMHigh)-calibratedEdges(250+FWHMLow);    
        FWHM(x,y) = FWHMHigh-FWHMLow; 

        peakCentroids(x,y) = peakCentroids(x,y)+250;
        if x==StepsX
            x=0;
            y=y+1;
        end    
    end

    fig = figure(12)
    imagesc((FWHM./peakCentroids)*100)
    colorMap = jet(256); 
    colormap(colorMap) 
    set(gca,'FontSize',18,'fontname','times') 
    caxis([0 3])
    d = colorbar; 
    d.Label.String = 'FWHM (%)'; 
    d.Label.FontSize = 18; 
    xlabel('X (Pix)','FontSize', 18) 
    ylabel('Y (Pix)','FontSize', 18) 
    axis square 
%     saveas(fig,'FWHMPercent')
    
    save('saveFWHMkeV','FWHM')
    FWHMPercent = ((FWHM./peakCentroids)*100);
    save('saveFWHMpercent','FWHMPercent')
end
