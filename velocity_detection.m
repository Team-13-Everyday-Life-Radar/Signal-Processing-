chirps = 16; % chirps per frame
samples = 128; % samples per chirp
receivers = 1; % text write and read does not yet work for receivers = 2
frames = 100; % how many frames do you want recorded?
current_line = 0;
j = 1; % frame counter in loops
h = 1; % chirp counter in loops
g = 1; % chirp counter in loops
PhDiff_array = zeros(1,chirps/2);
vr_array = zeros(1,chirps/2);
PhDiff_avg_array = zeros(1,frames);
vr_avg_array = zeros(1,frames);
% PhDiff_arrayRx2 = zeros(1,chirps/2);
% vr_arrayRx2 = zeros(1,chirps/2);
% PhDiff_avg_arrayRx2 = zeros(1,frames);
% vr_avg_arrayRx2 = zeros(1,frames);

% % plotting
sample_timeVec = 1:1:samples;
x1 = zeros(samples,1,chirps/2,frames);
x2 = zeros(samples,1,chirps/2,frames);

% Radial Velocity\T1256_16_slowFWD_foil - 30 frames

store_PhDiff = zeros(1,frames);
store_vr = zeros(1,frames);

nfft = 128;
c = 3e8;
fs = 640000;
ts = (300e-6)/256;
t =  0 : 1/fs : (ts)-(1/fs); 
f2 = 24.025e9;
f1 = 24.225e9;
chirpslope = (200e6)/(300e-6);
minrange = 0.95; maxrange= 3;
minbeatfreq_bleh = ((chirpslope)*2*minrange)/c 
maxbeatfreq_bleh = ((chirpslope)*2*maxrange)/c 
% minbeatfreq = (2*minrange*(200e6))/((3e8)*(300e-6));
% maxbeatfreq = (2*maxrange*(200e6))/((3e8)*(300e-6));
freq = (fs/nfft)*(-nfft/2:(nfft/2)-1);
fr = freq((nfft/2 +1):(end));

% % % Read ydata matrix from text file (can only upload data 1 receiver at a time)
% % Read Rx1 & Rx2 data in ydataRx1.txt and ydataRx2.txt
ydata_readRx1 = zeros(samples,receivers,chirps,frames);
% ydata_readRx2 = zeros(samples,receivers,chirps,frames);
input1 = load('zeroRx2_128_16.txt');
% input2 = load('Car\T1_128_16_20.txt');
for j = 1:frames
    for h = 1:(chirps)
        for g = 1:samples
            ydata_readRx1(g,1,h,j) = input1(g+current_line,1) + 1i*input1(g+current_line,2);
    %         ydata_readRx2(g,1,h,j) = input2(g+current_line,1) + 1i*input2(g+current_line,2);
        end
        current_line = current_line+g; 
    end
end
j = 1;

while j <= frames



for h = 1:2:chirps-1
    x1(:,1,f,j) = sqrt(((imag(ydata_readRx1(:,1,h,j))).^2) + ((real(ydata_readRx1(:,1,h,j))).^2));
    x2(:,1,f,j) = sqrt((imag(ydata_readRx1(:,1,h+1,j)).^2) + ((real(ydata_readRx1(:,1,h+1,j))).^2));
%     x1(:,1,f,j) = ydata_readRx1(:,1,h,j);
%     x2(:,1,f,j) = ydata_readRx1(:,1,h+1,j);
%     
    x1(:,1,f,j) = x1(:,1,f,j) - mean(x1(:,1,f,j));
    x2(:,1,f,j) = x2(:,1,f,j) - mean(x2(:,1,f,j));
    
    X1 = fftshift(fft(x1(:,1,f,j)));
    X2 = fftshift(fft(x2(:,1,f,j)));
    
    X1_half = X1(((nfft/2)+1):end);
    X2_half = X2(((nfft/2)+1):end);
    
    stepfreq = fr(2) - fr(1); % affected by fs and nfft
    minindx = round((minbeatfreq - fr(1))/stepfreq);
    maxindx = round((maxbeatfreq - fr(1))/stepfreq);
    

    
    X1_final = X1_half((minindx+1):(maxindx+1));
    X2_final = X2_half((minindx+1):(maxindx+1));
    
    [pk1, indx1] = max((X1_final));
    [pk2, indx2] = max((X2_final));
    PhDiff_array(1,a) = angle(X2(indx2)) - angle(X1(indx1)); % in radians
    vr_array(1,a) = (PhDiff_array(1,a)*((3e8)/((f2+f1)/2)))/(4*pi*300*(10^-6));
    f = f+1;
    a = a+1;
end 


    disp(j)
    j = j+1;  
    PhDiff_avg_array(1,j) = mean(PhDiff_array);
    vr_avg_array(1,j) = mean(vr_array)*3.6;
    disp(PhDiff_avg_array(1,j))
    disp(vr_avg_array(1,j))

end

max_val = max(vr_avg_array)
min_val = min(vr_avg_array)
avg_val = mean(vr_avg_array)
% x_ph = 0:(500e-6)*chirps:((500e-6)*(frames*chirps))-((500e-6*frames*chirps)*(1/frames))
% x_ph = 0:(1/frames)*4:4-((1/frames)*4)
x_ph = 0:1:frames;

plot(x_ph,vr_avg_array)
title('Radial Velocity')
xlabel("frames")
ylabel("velocity (km/h)")

plot(x_ph,PhDiff_avg_array)
title("Avg Phase Difference")
xlabel("frames")
ylabel("phase (radians)")



