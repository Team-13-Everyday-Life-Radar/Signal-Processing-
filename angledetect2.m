clc
disp('******************************************************************');
addpath('..\..\RadarSystemImplementation'); % add Matlab API
clear all %#ok<CLSCR>
close all
resetRS; % close and delete ports

antenna_spacing = 6.22e-3; % in meters
c0 = 3e8; % Speed of light in vacuum

% 1. Create radar system object
szPort = findRSPort; % scan all available ports
oRS = RadarSystem(szPort); % setup object and connect to board

disp('Connected RadarSystem:');
oRS %#ok<*NOPTS>

% 2. Enable automatic trigger with frame time 1s
% oRS.oEPRadarBase.set_automatic_frame_trigger(1000000);

 oRS.oEPRadarBase.stop_automatic_frame_trigger; % stop it to change values
%
%    Properties of radar
%
   oRS.oEPRadarFMCW.lower_frequency_kHz = 24025000; % lower FMCW frequency 
   oRS.oEPRadarFMCW.upper_frequency_kHz = 24225000; % upper FMCW frequency
   oRS.oEPRadarFMCW.tx_power = 4;
   oEPTargetDetection.max_range_cm = 200;
   oRS.oEPRadarBase.num_chirps_per_frame = 16;  %[1,128]
   oRS.oEPRadarBase.num_samples_per_chirp = 64; % [32, 64, 128, 256] 
   oRS.oEPRadarADCXMC.samplerate_Hz = 640000; % 640000 affects the chirp time
   oRS.oEPRadarFMCW.direction = 'Up Only';
   oRS.oEPRadarBase.set_automatic_frame_trigger(1000000);
   
   chirp = 16;
   i = 1;
   n = 1;
   angle_sum = 0;
   
while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data
    
%    disp(ydata);
    
    %set value for lambda
    fC = ((24225000 + 24025000)/2)*1e3;
    lambda = c0/fC;
    angle_sum = 0;
       
        for i = 1:chirp
     
            ang_Rx1 = angle(ydata(:,1,i)); % Phase of received signal for Rx1
            ang_Rx2 = angle(ydata(:,2,i)); % Phase of received signal for Rx2
                      
            d_phi = mean(ang_Rx1) - mean(ang_Rx2); % Phase difference between the received signal at the two receivers 
           
            target_angle = (asin((d_phi * lambda) / (antenna_spacing * (2*pi)))); % AOA in radians
    
            target_angle_deg = ((target_angle) * 180 / pi); % AOA in degrees
            angle_sum = angle_sum + target_angle_deg;

        end
        angle_avg = angle_sum / chirp

end