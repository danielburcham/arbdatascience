function [C,fig_handle] = parsedx(filename)
%{
PARSEDX
---------------------------------------------------------------------------
DESCRIPTION
Parsedx accepts an Excel file containing strings recorded from a Dillon 
EDXtreme dynamometer and returns data parsed into separate datetime, 
elapsed seconds, instantaneous force, and peak force columns. The strings
are printed by the Dillon EDXtreme in continuous recording mode in #4 
format with the dynamometer connected to a computer by a serial port 
connection and WedgeLink software. 

INPUTS
filename: string - filepath of Excel file containing original strings

OUTPUTS
B: matrix - m x 4 matrix containing datetime, elapsed seconds,
instantaneous force, and peak force records in separate columns
fig_handle: figure object - line plot of elapsed seconds and instantaneous
force
---------------------------------------------------------------------------
%}
% Read file
[~,A,~]=xlsread(filename);
% Initialize array
B = zeros(size(A,1),2);
% Read and convert datetime stamps
dt = regexp(A,'\d*\s\w*\s\d{4}[\,]\d*\:\d{2}\:\d{2}','match');
t = datetime(string(dt),'InputFormat','d MMM y,H:mm:ss');
s = seconds(t-t(1));
% Read and convert force measurements
f = regexp(A,'\d*\.\d{1}','match');
for i = 1:size(f,1)
    B(i,1:2) = [str2double(f{i,1}{1,1}) str2double(f{i,1}{1,2})];
end
C = table(t,s,B(:,1),B(:,2),'VariableNames',{'Datetime','Seconds',...
    'F_ins','F_pk'});
% Plot
fig_handle = plot(s,B(:,1));
xlabel('Time, t (sec)');
ylabel('Force, F (kg)');
writetable(C,filename,'sheet',2);