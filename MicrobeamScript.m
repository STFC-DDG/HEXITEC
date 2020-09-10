clear
clc
tic

addpath C:\Users\dml36958\Desktop\DesktopDocs\MatlabFuncs %Any path where read80x80 function saved!

Home_Path = 'F:\';
cd(Home_Path)
Home_Folder = cd;

directories_list = dir('3348*.*');

goodRegionScans = [66 67 68 69 70 71 72 74]';

for i=1:length(goodRegionScans)
    for j = 1:length(directories_list)
        name = strcat('3348', num2str(goodRegionScans(i)));
        strs ={directories_list.name};
        index = find(ismember(strs, name));
        goodDirectoriesList(i) = directories_list(index);
    end
end


%%
numberSkipped = 0;
StepsY = 41;
StepsX = 41;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***PIXEL POSITION***
pxPosXmin = 14;
pxPosXmax = 21;
pxPosYmin = 9;
pxPosYmax = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binWidth = 1;
binMax = 8000;
bins = (0:binWidth:binMax)';
spectrasRaw=zeros(StepsX,StepsY,binMax);
spectrasCSD=zeros(StepsX,StepsY,binMax);


pxArray = zeros(StepsX,StepsY);
pxArrayHxt = zeros(StepsX,StepsY);
storedFrames = zeros(StepsX,StepsY);
normPxArray = zeros(StepsX, StepsY);
discriminatedPxArray = zeros(StepsX, StepsY);
disArray = zeros(StepsX, StepsY);
normDisArray = zeros(StepsX, StepsY);
shareIm = zeros(StepsX, StepsY);

fakeHxt = zeros(8,7,binMax*(1/binWidth),(StepsY*StepsX));
shareData = zeros(8,7,binMax*(1/binWidth),(StepsY*StepsX));
disData = zeros(8,7,binMax*(1/binWidth),(StepsY*StepsX));
numShareFrames = 0;
timePred=0;

j=1;
x=1;

for CD=1:length(goodDirectoriesList)
    cd (strcat(goodDirectoriesList(CD).folder,goodDirectoriesList(CD).name))
    currentSubDirs = dir('*.*');
    
    currentSubDirs = currentSubDirs(3:length(currentSubDirs));


%     cd G:\334890    
%     CdArraysListHxt = dir('**/*.hxt');
       
    
    %sets number of values to read in so as to remove any overlapping data
    %i.e. reads only complete columns of data
    loopCutOff = length(currentSubDirs);    
    if length(currentSubDirs)>StepsX
        loopCutOff = (floor(length(currentSubDirs)/StepsY))*StepsY;
    end
    
    i=1;    
    for count=1:loopCutOff
        cd (strcat(currentSubDirs(count).folder,'\',currentSubDirs(count).name))
        CdArraysList = dir('**/*.bin');        
        fullPath = strcat(CdArraysList.folder, '\', CdArraysList.name);
        
        testFull = size(CdArraysList);
        if testFull(1) == 0
            continue
        end
        FileSize = CdArraysList.bytes;        
        
        fid = fopen(fullPath);
        fseek(fid,0,'bof');
        nFrames = floor(FileSize/(80*80));

        currentSpectra=zeros(length(bins),1);
        currentSpectraCSD=zeros(length(bins),1);
        
        for fNum = 1:((nFrames/2)-1)       
            frame = read80x80Frame(fid,1);            
            frameMask = frame>100;
            frame = immultiply(frame,frameMask);

            binnedFrame = discretize(frame(pxPosXmin:pxPosXmax,pxPosYmin:pxPosYmax),bins);
            binnedFrame=binnedFrame(binnedFrame>1);
            for z=1:length(binnedFrame)
                currentSpectra(binnedFrame(z))=currentSpectra(binnedFrame(z))+1;
            end                
                        
            clusters = bwlabel(frameMask(pxPosXmin:pxPosXmax,pxPosYmin:pxPosYmax));
            clusterMask = clusters>0;
            numClusters = max(max(clusters));
            frameCSD = immultiply(frame(pxPosXmin:pxPosXmax,pxPosYmin:pxPosYmax),clusterMask);
            for currentCluster=1:numClusters
                [x,y]=find(clusters==currentCluster);
                if length(x)>1
                    if currentCluster ==1
                        numShareFrames = numShareFrames + 1;
                    end
                    for z=1:length(x)
                        frameCSD(x(z),y(z))=0;
                    end
                end
            end
            binnedFrameCSD = discretize(frameCSD,bins);
            binnedFrameCSD=binnedFrameCSD(binnedFrameCSD>0);
            for z=1:length(binnedFrameCSD)
                currentSpectraCSD(binnedFrameCSD(z))=currentSpectraCSD(binnedFrameCSD(z))+1;
            end       
        end
        spectrasCSD(i,j,:)=squeeze(currentSpectraCSD(2:8001)); 
        spectrasRaw(i,j,:)=squeeze(currentSpectra(2:8001));
        
        fclose(fid);      
        storedFrames(i,j) = (nFrames/2)-1;     
              
        timePred = timePred+1;
        if timePred == 10
            time10 = toc;
            timePrediction = (time10/10)*(StepsX*StepsY);
            totalMinutes = timePrediction/60
        end
        
        i = i+1;

        if i >StepsY
            j = j+1
            i = 1;
        end    
        if j >StepsX-1
            j=1
        end
    end
    cd (Home_Path)      
end


%%
cd ('C:\Users\dml36958\Desktop\BeamTime_FinalRun\Bad12keV')
ImRaw = sum(spectrasRaw,3);
fig = figure(1)
ImToPlot = (ImRaw./storedFrames)/max(max((ImRaw./storedFrames)));
imagesc(ImToPlot(:,1:22))
% colorMap = jet(256); 
% colormap(colorMap) 
set(gca,'FontSize',18,'fontname','times') 
d = colorbar; 
d.Label.String = 'Normalised Events'; 
d.Label.FontSize = 18; 
xlabel('X (Pix)','FontSize', 18) 
ylabel('Y (Pix)','FontSize', 18) 
axis square 
caxis([0.2 0.95])
saveas(fig,'RawIm.fig')

ImCSD = sum(spectrasCSD,3);
fig = figure(2)
normToFrames = (ImCSD/max(max((ImCSD))));
% imagesc((ImCSD./storedFrames)/max(max((ImRaw./storedFrames))))
imagesc(normToFrames(:,1:40))
caxis([0.1 1])
% colorMap = jet(256); 
% colormap(colorMap) 
set(gca,'FontSize',18,'fontname','times') 
d = colorbar; 
d.Label.String = 'Normalised Events'; 
d.Label.FontSize = 18; 
xlabel('X (Pix)','FontSize', 18) 
ylabel('Y (Pix)','FontSize', 18) 
axis square 
% saveas(fig,'CSDIm.fig')

%%

y=1;
x=0;
peakAmplitudes=zeros(StepsX,StepsY);
peakCentroids=zeros(StepsX,StepsY);
peakHalf=zeros(StepsX,StepsY);
FWHM=zeros(StepsX,StepsY);

for z = 1:(StepsX*StepsX)
    x=x+1;

    currentSubPix = squeeze(spectrasCSD(x,y,:));
          
    peakAmplitudes(x,y)=max(currentSubPix);
    peakCentroids(x,y)=find(currentSubPix==peakAmplitudes(x,y),1,'first');
    peakHalf(x,y) =  peakAmplitudes(x,y)/2;
    FWHMHigh = peakCentroids(x,y) + (find(currentSubPix(peakCentroids(x,y):length(currentSubPix))< peakHalf(x,y),1,'first'));
    FWHMLow = find(currentSubPix>peakHalf(x,y), 1, 'first');
    if isempty(FWHMHigh)==1
        continue
    end
    FWHM(x,y) = FWHMHigh-FWHMLow;    
 
    if x==StepsX
        x=0;
        y=y+1;
    end
end

fig = figure(3)
imagesc(peakCentroids(:,1:40))
% colorMap = jet(256); 
% colormap(colorMap) 
set(gca,'FontSize',18,'fontname','times') 
d = colorbar; 
d.Label.String = 'Energy (ADU)'; 
d.Label.FontSize = 18; 
xlabel('X (Pix)','FontSize', 18) 
ylabel('Y (Pix)','FontSize', 18) 
axis square 
caxis([100 350])
% saveas(fig,'EnergyMapADU.fig')

fig = figure(4)
imagesc(FWHM(:,1:40))
% colorMap = jet(256); 
% colormap(colorMap) 
set(gca,'FontSize',18,'fontname','times') 
d = colorbar; 
d.Label.String = 'FWHM (ADU)'; 
d.Label.FontSize = 18; 
xlabel('X (Pix)','FontSize', 18) 
ylabel('Y (Pix)','FontSize', 18) 
axis square 
caxis([30 150])
% saveas(fig,'FWHM.fig')

% save('matlab')
