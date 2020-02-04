
clc
disp('******************************************************************');
addpath('C:\IFX_P2G-HW-SW_V1.0.2\Firmware_Software\Communication Library\ComLib_Matlab_Interface\RadarSystemImplementation'); % add Matlab API
clear all; %#ok<CLSCR>
close all
resetRS; % close and delete ports

% 1. Create radar system object
szPort = findRSPort;  %scan all available ports
oRS = RadarSystem(szPort); % setup object and connect to board

disp('Connected RadarSystem:');
oRS %#ok<*NOPTS>

% 2. Enable automatic trigger with frame time 1s
oRS.oEPRadarBase.set_automatic_frame_trigger(1000000);

while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data
    
     disp((ydata));
     
    tf =  double(oRS.oEPRadarBase.chirp_duration_ns);% micro 
    t = linspace(0,1e-3.*tf,64); %frame + 1; %linspace(0,64,64);
    fc = double(oRS.oEPRadarDoppler.frequency_kHz);


    A= sqrt(real(ydata(:,1)).^2 + imag(ydata(:,1)).^2);    
    angle= atan(imag(ydata(:,1))./ real(ydata(:,1)));
    r = real(ydata(:,1));
    j = imag(ydata(:,1));
    I = r'.*cos(2.*pi*fc*t);
    Q = j'.*sin(2.*pi*fc*t);

figure(1)
    plot(t, I,'b')
    hold on 
    plot(t,Q,'r')
    grid on
    title('I/Q Signal');
    xlabel('time (ns)');
    ylabel ('amplitude');
    legend('I signal', 'Q signal');
   
 figure(2)
    plot(t,A'.*cos(2.*pi.*fc.*t + angle'));% - Q.*sin(2.*pi.*10000000.*frame))
    grid on
    title('IF Signal/ Receiver Output');
    xlabel('time (ns)');
    ylabel ( 'amplitude');


end
