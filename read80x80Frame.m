function [frame] = read80x80Frame(fid,flag)
%Function to return 3 frames of 80x80 together in a 80x240 frame of data
%for analysis of HEXITEC .bin file data.

%fid is the file identifier generated from fopen in main program
%flag to in indicate if the data should be sorted into the correct physical
%orientation or with the raw interleaved format for speed.

% fid = fopen('C:\TestData8\Test.bin');
% flag = 1;

%Matt Wilson 14/07/16
thisID = fid;
frame = zeros(80,80);  %create empty frame
if ~feof(fid)
    if flag == 1
        %Sort Data into 240x80 pixels with correct orientation
        dataIn = fread(fid,6400,'uint16');
        n = 1;
        for row=1:80
            for j=1:20
                for k=1:4
                    col = j+((k-1)*20);
                    pixVal = dataIn(n);
                    frame(row,col) = pixVal(1,1);
                    n= n+1;
                end
            end
        end
        
        
        
    elseif flag == 0
        %Read Data into 240x80 pixels with interleaved format for speed
        frame = reshape(fread(fid,6400,'uint16'),80,80);
        
    else
        disp('Invalid Flag. Use Flag = 1 for sorting frame into correct orientation or Flag = 0 for interleaved frame which will be faster.')
    end
else
    disp('End of File')
end

% end