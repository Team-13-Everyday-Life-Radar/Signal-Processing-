%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = extract_raw_data (in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (c) 2014-2019, Infineon Technologies AG
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,are permitted provided that the
% following conditions are met:
%
% Redistributions of source code must retain the above copyright notice, this list of conditions and the following
% disclaimer.
%
% Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
% disclaimer in the documentation and/or other materials provided with the distribution.
%
% Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote
% products derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE  FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
% WHETHER IN CONTRACT, STRICT LIABILITY,OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% This simple example demos the acquisition of data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cleanup and init
% Before starting any kind of device the workspace must be cleared and the
% MATLAB Interface must be included into the code. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
disp('******************************************************************');
addpath('..\..\RadarSystemImplementation'); % add Matlab API
clear all %#ok<CLSCR>
close all
resetRS; % close and delete ports

% 1. Create radar system object
szPort = findRSPort; % scan all available ports
oRS = RadarSystem(szPort); % setup object and connect to board

disp('Connected RadarSystem:');
oRS %#ok<*NOPTS>

% 2. Enable automatic trigger with frame time 1s
oRS.oEPRadarBase.set_automatic_frame_trigger(1000000);
% 
% oRS.oEPRadarBase.stop_automatic_frame_trigger; % stop it to change values 
%  
%   oRS.oEPRadarFMCW.lower_frequency_kHz = 24025000; % lower FMCW frequency   
%   oRS.oEPRadarFMCW.upper_frequency_kHz = 24225000; % upper FMCW frequency   
%   oRS.oEPRadarFMCW.tx_power = oRS.oEPRadarBase.max_tx_power;   
%   oRS.oEPRadarBase.num_chirps_per_frame = 1;   
%   oRS.oEPRadarBase.num_samples_per_chirp = 256; % [32, 64, 128, 256]   
%   oRS.oEPRadarBase.rx_mask = bin2dec('0011'); % enable two RX antennas   
%   oRS.oEPRadarFMCW.direction = 'Up Only';

% Tc = cast(oRS.oEPRadarBase.chirp_duration_ns * (10^-9),'double'); % chirp time
oRS.oEPTargetDetection.min_range_cm = 300; % set max distance
oRS.oEPTargetDetection.max_speed_kmh = 10; % set max speed
i = 0;

% % Variables needed to write ydata to text file
% Change variables accordingly
chirps = 16; % How many chirps per frame?
samples = 64; % How many samples per chirp?
frames = 5; % How many consecutive frames do you want recorded?
receivers = 1; % fixed variable,i.e., don't change
current_line = 0; % fixed variable
j = 1; % fixed variable
fileID1 = fopen('ydataRx1.txt','w'); % Receiver 1 ydata
fileID2 = fopen('ydataRx2.txt','w'); % Receiver 2 ydata

while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data

    % % % % Write ydata matrix to text file (can only upload data from 1 receiver at a time)
    % % Write data from Rx1
    if j <= frames
        for i = 1:chirps
            fprintf(fileID1,'%f %f\n',[real(ydata(:,1,i,j)),imag(ydata(:,1,i,j))].'); 
        end
    end
    % % Write data from Rx2
    if j <= frames
        for i = 1:chirps
            fprintf(fileID2,'%f %f\n',[real(ydata(:,2,i,j)),imag(ydata(:,2,i,j))].'); 
        end
    else
        break
    end
    j = j+1;
    
end

% % Conclude write to text file for receiver 1 and 2
fclose(fileID1); 
fclose(fileID2); 
