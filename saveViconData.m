clear
close all

dbstop if error

addpath Vicon_Data





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg stop at corner clockwise
test1.VICON_FILENAME = 'Subj05_5deg_Test1.csv';
test1.testObj = getVicondata(test1.VICON_FILENAME);
test1.testNum = 1;
rangeVicon = 1:23240;
test1.rangeVicon = rangeVicon;
disp('finish: 1')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg stop at corner counterclockwise
test2.viconFilename = 'Subj05_5deg_Test2.csv';
test2.testObj = getVicondata(test2.viconFilename);
test2.testNum = 2;
rangeVicon = 1:23840;
% rangeVicon = 5132:1:5500;
test2.rangeVicon = rangeVicon;
disp('finish: 2')

% testDataAll = {test2};
% disp('Saving...')
% save('testDataAll.mat','testDataAll')
% disp('Save done.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg non stop clockwise
test3.viconFilename = 'Subj05_10deg_Test3.csv';
test3.testObj = getVicondata(test3.viconFilename);
test3.testNum = 3;
% rangeVicon = 5054:1:5310;
rangeVicon = 1:22920;
test3.rangeVicon = rangeVicon;
disp('finish: 3')

% testDataAll = {test3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg non stop counterclockwise
test4.viconFilename = 'Subj05_10deg_Test4.csv';
test4.testObj = getVicondata(test4.viconFilename);
test4.testNum = 4;
rangeVicon = 1:23710;
% rangeVicon = 18829:1:18976; % one step
test4.rangeVicon = rangeVicon;
disp('finish: 4')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20 deg stop at corner clockwise
test5.viconFilename = 'Subj05_15deg_Test5.csv';
test5.testObj = getVicondata(test5.viconFilename);
test5.testNum = 5;
rangeVicon = 1:23450;
% rangeVicon = 20200:1:20345; % one step
test5.rangeVicon = rangeVicon;
disp('finish: 5')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20 deg stop at corner counterclockwise
test6.viconFilename = 'Subj05_15deg_Test6.csv';
test6.testObj = getVicondata(test6.viconFilename);
test6.testNum = 6;
rangeVicon = 1:23400;
% rangeVicon = 21195 :1:21323; % one step
test6.rangeVicon = rangeVicon;
disp('finish: 6')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20 deg non stop clockwise
test7.viconFilename = 'Subj05_20deg_Test7.csv';
test7.testObj = getVicondata(test7.viconFilename);
test7.testNum = 7;
rangeVicon = 1:27200;
% rangeVicon = 20200:1:20345; % one step
test7.rangeVicon = rangeVicon;
disp('finish: 7')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20 deg non stop counterclockwise
test8.viconFilename = 'Subj05_20deg_Test8.csv';
test8.testObj = getVicondata(test8.viconFilename);
test8.testNum = 8;
rangeVicon = 1:23780;
% rangeVicon = 21195 :1:21323; % one step
test8.rangeVicon = rangeVicon;
disp('finish: 8')



testDataAll = {test1, test2, test3,test4, test5, test6, test7, test8};

% testDataAll = {test1, test2, test3, test4};

disp('Saving...')
save('testDataAll.mat','testDataAll')
disp('Save done.')

