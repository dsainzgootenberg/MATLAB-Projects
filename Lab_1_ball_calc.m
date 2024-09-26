clear
close all
clc

% Sets up a the array showing the location of where to read images from   
% with the folder var allowing for easy switching to the different
% collec tions of images%
Folder = {'\Pickle Ball';'\Ping Pong Ball';'\Squash Ball';'\Tennis Ball'};
%Framerate
fr = 2000;
%pixels to in. conversion coeff
px2in = 37;
for m = 1:length(Folder)
    theFiles = dir(fullfile('C:\Users\Dany SG\Documents\CWRU Classes\Fall 2024\EMAE 285\Labs\Full reports\MATLAB', Folder{m}, '*.tif'));

    %Creates a Structure that houses all filemanes and Image arrays as well as
    %any potential future calcs that will be needed
    image_Data = struct('filename', {}, 'Image', {}, 'X', {}, 'Y', {});

    %Creates a for loop to read all the images in the specified folder and
    %place any relevant data into the imageData structure
    for k = 1:10:length(theFiles)
        k
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(theFiles(k).folder, baseFileName);
        Image = imread(fullFileName);
        image_Data(k).filename = baseFileName;
        if k == 1
            image(Image)
            % Makes user draw line to find potential diameter of the ball and calulates
            % the radius within +-5% to account for drawing accuracies
            d = drawline;
            pos = d.Position;
            diffPos = diff(pos);
            Dia = hypot(diffPos(1), diffPos(2));
            radMin = round(Dia/2 * 0.90);
            radMax = round(Dia/2 * 1.1);
        end
        % Detects and marks a circle using max an min radi calculated to ensure only one is
        % found and returns its center and radi where the center will be collected
        % and stored with center of other images
        [ c , r ] = imfindcircles(Image , [radMin , radMax] , Sensitivity=0.95 , ObjectPolarity="dark" , EdgeThreshold=0.1 );
        while isempty(c)
            [ c , r ] = imfindcircles(Image , [radMin , radMax] , Sensitivity=0.98 , ObjectPolarity="dark" , EdgeThreshold=0.1 );
            if isempty(c)
                image(Image)
                d = drawline;
                pos = d.Position;
                diffPos = diff(pos);
                Dia = hypot(diffPos(1), diffPos(2));
                radMin = round(Dia/2 * 0.90);
                radMax = round(Dia/2 * 1.1);
            end
        end
        image_Data(k).X = c(1 , 1);
        image_Data(k).Y = c(1 , 2);

        %plots a 5 x 5 dot in the center for tracking
        Image( round(c(2))-2:round(c(2))+2, round(c(1))-2:round(c(1))+2 ) = 255;
        image(Image)
        radMin = round(r(1) * 0.9);
        radMax = round(r(1) * 1.1);
        cir = viscircles(c, r);
        image_Data(k).Image = Image;
    end
    h1(m) = {single([image_Data.Y])'};
    
end
save('H.mat', 'h1');
hold off
for l=1:4
t = [1:length(h1{l})]./fr
h_in = [1024-h1{l}]./37
v(l) = plot(t, h_in, 'b-', 'linewidth', 3)
hold on
end
hold off
title('Height of balls per Time')
xlabel('Time (s)') 
ylabel('Hieght (in)')
legend({'Pickel Ball','Ping Pong Ball','Squash Ball','Tennis Ball'},'Location','Northeast')
    shg
