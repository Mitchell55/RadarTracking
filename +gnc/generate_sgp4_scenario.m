function [sc, sat, sat_ID] = generate_sgp4_scenario(startTime,stopTime,dt,tleFile,OP)
% This function creates a satellite scenario and adds satellites from tle
% file and ground stations.
%
% Inputs: startTime - simulation start time (datetime format)
%         stopTime - simulation end time (datetime format)
%         dt - step size (sec)
%         tleFile - TLE file name
%         OP - orbit propagator string ('SGP4')
%
% Outputs: sc - satelliteScenario object (read only)
%          sat - satellite object(read only)
%          sat_ID - satellite ID from TLE (index cooresponds to sat ID
%                   assigned by satellite object
%          gs - ground station object
% 
% Mikaela Dobbin 06/21/2021
% Laura Davies 2/07/2022

% Create satellite scenario object
sc = satelliteScenario(startTime, stopTime, dt);

% Add satellites to satellite scenario using sgp4 propagator
sat = satellite(sc, tleFile, 'OrbitPropagator', OP);

% Save sat ID from tle for output file
inFile = fopen(tleFile, 'r'); 
line = true;
count = 1;
while (~feof(inFile))
    longstr = fgets(inFile);
    if length(longstr)>30 && line==true
        ID = string(longstr(3:7));
        sat_ID(count,1) = str2double(ID);
        line = false;
        count = count + 1;
    elseif length(longstr)>30 && line == false
        line = true;
    end       
end
fclose(inFile);
end