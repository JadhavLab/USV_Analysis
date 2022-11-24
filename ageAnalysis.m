% ageAnalysis
%{
this analysis is basically to examine some surface level characteristics of
the calls and plot them out.  We know that the number of calls, their
intensity and duration all trail off as the animals grow up.  We also know
that they increase in complexity.  So the goal here is to examine this
transition- which parameters change and when- do they all change at once or
do they transition separately?

The other surface level question is are these parameters different either
early or late across sex or genotype, and if not, does this transition
occurr differently for different sexes or genotypes
%}



%% first lets gather some call parameters
callFolder='E:\Brandeis datasets\FMR1 Project Data\USV data\segData';
for i=1:length(USVlegend)
    callTimes=USVlegend.callTimes(i);
    dataFile= load(fullfile(callFolder,USVlegend.name(i)),'blobs','segCalls','params','spect');


end