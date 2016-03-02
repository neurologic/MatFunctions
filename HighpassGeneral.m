function [data]=HighpassGeneral(data,fs)
%data must be in format rows=trials columns=data
% data=BaselineGeneral(data,baselinewin,basedata);
% shouldn't need to baseline because high pass removes DC offset

%cutoff frequency hard-coded at 100 here cause i have never used anything
%different ever
data=fftFilter(data',fs,100,2)';
