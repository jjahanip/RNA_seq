function setPath()
% setPath()
%
% Update the search path with all the directories
% necessary for running the algorithm.
%
% 

basepath = cd();
addpath(genpath([basepath, '\lib']));



