code = [0 1 0 1 0 1 0 1 0 1;
    0 1 0 1 0 1 0 1 0 1;
    0 1 0 1 0 1 0 1 0 1;
    0 1 0 1 0 1 0 1 0 1;
    0 1 0 1 0 1 0 1 0 1;
    0 1 0 1 0 1 0 1 0 1;
        %1 0 1 0 1 0 1 0 1 0;
        0 1 0 1 0 1 0 1 0 1;
        %1 0 1 0 1 0 1 0 1 0;
        0 1 0 1 0 1 0 1 0 1;
        %1 0 1 0 1 0 1 0 1 0;
        0 1 0 1 0 1 0 1 0 1;
        %1 0 1 0 1 0 1 0 1 0;
        0 1 0 1 0 1 0 1 0 1;];
        %1 0 1 0 1 0 1 0 1 0;];
t = 3;
% MATLAB Code from Sensor Array Analyzer App

% Generated by MATLAB 9.9 and Phased Array System Toolbox 4.4

% Generated on 16-Apr-2021 20:03:39

% Create a uniform rectangular array
[row, column] = size(code);
Array = phased.URA('Size',[row column],...
'Lattice','Rectangular','ArrayNormal','x');
% The multiplication factor for lambda units to meter conversion
Array.ElementSpacing = [1/t 1/t]*0.001;
% Calculate Row taper
rwind = ones(1,row).';
% Calculate Column taper
cwind = ones(1,column).';
% Calculate taper
taper = rwind*cwind.';
Array.Taper = taper.';

% Create an isotropic antenna element
Elem = phased.IsotropicAntennaElement;
Elem.FrequencyRange = [0 300000000000];
Array.Element = Elem;

% Partition the array
Array = phased.PartitionedArray('Array',Array,...
 'SubarraySelection',rewriteCode(code));
Array.SubarraySteering = 'Phase';
Array.NumPhaseShifterBits = 1000;
Array.PhaseShifterFrequency = 300000000000;
% Assign Frequencies and Propagation Speed
Frequency = 300000000000;
PropagationSpeed = 300000000;

SubarraySteerAngles = [0;0];

% Create Figure

% Plot Array Geometry
figure;
viewArray(Array,'ShowNormal',true,...
  'ShowTaper',false,'ShowIndex','None',...
  'ShowSubarray','All');

% Calculate Steering Weights

Freq3D = 300000000000;
% Find the weights
w = ones(getNumSubarrays(Array), length(Frequency));

% Plot 3d graph
format = 'polar';
figure;
pattern(Array, Freq3D , 'PropagationSpeed', PropagationSpeed,...
 'SteerAngle', SubarraySteerAngles, ...
 'Type','directivity', 'CoordinateSystem', format,'weights', w(:,1));

% Find the weights
w = ones(getNumSubarrays(Array), length(Frequency));

% Plot 2d azimuth graph
format = 'polar';
cutAngle = 0;
figure;
pattern(Array, Frequency, -180:180, cutAngle, 'PropagationSpeed', PropagationSpeed,...
 'SteerAngle', SubarraySteerAngles, ...
 'Type', 'directivity', 'CoordinateSystem', format ,'weights', w);

% Find the weights
w = ones(getNumSubarrays(Array), length(Frequency));

% Plot 2d elevation graph
format = 'polar';
cutAngle = 0;
figure;
pattern(Array, Frequency, cutAngle, -90:90, 'PropagationSpeed', PropagationSpeed,...
 'SteerAngle', SubarraySteerAngles, ...
 'Type', 'directivity', 'CoordinateSystem', format ,'weights', w);

function array = rewriteCode(code)
    [r, c] = size(code);
    array = zeros(2,r*c);
    for i = 1:r*c
        array(1,i) = -1*code(i) + 1;
        array(2,i) = code(i);        
    end
end

