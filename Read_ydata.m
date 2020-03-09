chirps = 16; % chirps per frame
samples = 64; % samples per chirp
receivers = 1; % text write and read doesn't work for receivers = 2
frames = 5; % how many frames do you want recorded?
current_line = 0;
j = 1;

% % Read ydata matrix from text file (can only upload data 1 receiver at a time)
% Read Rx1 & Rx2 data in ydataRx1.txt and ydataRx2.txt
ydata_readRx1 = zeros(samples,receivers,chirps,frames);
ydata_readRx2 = zeros(samples,receivers,chirps,frames);
input1 = load('ydataRx1.txt');
input2 = load('ydataRx2.txt');
while j <= frames
for h = 1:(chirps)
    for g = 1:samples
        ydata_readRx1(g,1,h,j) = input1(g+current_line,1) + 1i*input1(g+current_line,2);
        ydata_readRx2(g,1,h,j) = input2(g+current_line,1) + 1i*input2(g+current_line,2);
    end
    current_line = current_line+g;
end
j = j+1;
end
disp(ydata_readRx1) % Receiver 1 ydata stored in ydata_readRx1 array
disp(ydata_readRx2) % Receiver 2 ydata stored in ydata_readRx1 array
