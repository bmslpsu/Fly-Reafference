function [s,devices,ch] = daq_realtime(rate,AI,AO)
%% daq_realtime: creates DAQ session objects for MC_USB_1208FS_PLUS
%   INPUTS:
%       rate     	: aquisition rate
%       AI          : analog input channel #'s
%
%   OUTPUTS:
%       s           : session objects
%       devices    	: DAQ devices object
%       ch          : channels
%

devices = daq.getDevices;
ID = {devices.ID};
readID = ID{1};
s.Read = daq.createSession('mcc'); % create session to read inputs

s.Read.Rate = rate; % samples per second
s.Read.IsContinuous = true; % continuous data collection until stopped

% AI = [1 2 3]; % analog input channels
ch.AI = addAnalogInputChannel(s.Read, readID, AI, 'Voltage'); % add analog input channels

fprintf('MC_USB_1208FS_PLUS: \n Rate: %i \n AI: %s \n',...
    s.Read.Rate, num2str(AI))

end