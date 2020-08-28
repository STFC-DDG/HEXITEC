filesList = ['D212801'];

sqEdge = 80;
axis = 0:0.02:500;
globalHistCalibrated = zeros(1,length(axis));

for fNum=1:length(filesList(:,1))
    fNum

%     load(strcat('I:\Illinois\',filesList(fNum,:),'\Figures_',filesList(fNum,:),'\savepixelHistCSD.mat'))
%     load(strcat('I:\Illinois\',filesList(fNum,:),'\Figures_',filesList(fNum,:),'\saveCalibratedEdges.mat'))

    load(strcat('I:\Illinois\',filesList(fNum,:),'\Figures_Am241_D212801_300s_500V_200817_153836.bin','\savepixelHistCSD.mat'))
    load(strcat('I:\Illinois\',filesList(fNum,:),'\Figures_Am241_D212801_300s_500V_200817_153836.bin','\saveCalibratedEdges.mat'))

    for i =1:sqEdge
        for j=1:sqEdge
            for z=1:800
                currentBin = ceil(squeeze(calibratedEdges(z,i,j)) * 50) / 50;
                bin = find(axis==currentBin);

                if bin <= 0
                    continue
                end
                currentVal = squeeze(pixelHistCSD(z,i,j));
                globalHistCalibrated(bin) = globalHistCalibrated(bin) + currentVal;
            end
        end
    end

    nonZeros = find(globalHistCalibrated > 0);

    fig=figure
    plot(axis(nonZeros),smooth(globalHistCalibrated(globalHistCalibrated>0)),'b')   
    set(gca,'FontSize',18,'fontname','times') 
    xlabel('Energy (keV)','FontSize', 18)
    ylabel('Counts)','FontSize', 18)
    title('Global Spectra - CSD')
    box off
    saveas(fig,(strcat('I:\Illinois\',filesList(fNum,:),'\Figures_',filesList(fNum,:),'\CalibratedGlobalCSDplot.fig')))
    close all

    clear pixelHistCSD
    clear calibratedEdges
end