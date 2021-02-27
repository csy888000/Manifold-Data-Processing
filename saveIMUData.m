clear
close all

dbstop if error

addpath IMU_data





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg stop at corner clockwise
test1.IMU_FILENAME = '210224_221242_test_1.xls';
test1.IMUobj = getIMUdata(test1.IMU_FILENAME);
disp('finish: 1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg stop at corner counterclockwise
test2.IMU_FILENAME = '210224_222009_test_2.xls';
test2.IMUobj = getIMUdata(test2.IMU_FILENAME);
disp('finish: 2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg non stop clockwise
test3.IMU_FILENAME = '210224_222938_test_3.xls';
test3.IMUobj = getIMUdata(test3.IMU_FILENAME);
disp('finish: 3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg non stop counterclockwise
test4.IMU_FILENAME = '210224_223548_test_4.xls';
test4.IMUobj = getIMUdata(test4.IMU_FILENAME);
disp('finish: 4')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg stop at corner clockwise
test5.IMU_FILENAME = '210224_224721_test_5.xls';
test5.IMUobj = getIMUdata(test5.IMU_FILENAME);
disp('finish: 5')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg stop at corner counterclockwise
test6.IMU_FILENAME = '210224_225311_test_6.xls';
test6.IMUobj = getIMUdata(test6.IMU_FILENAME);
disp('finish: 6')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg non stop clockwise
test7.IMU_FILENAME = '210224_230117_test_7.xls';
test7.IMUobj = getIMUdata(test7.IMU_FILENAME);
disp('finish: 7')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 deg non stop counterclockwise
test8.IMU_FILENAME = '210224_230648_test_8.xls';
test8.IMUobj = getIMUdata(test8.IMU_FILENAME);
disp('finish: 8')


IMUDataAll = {test1, test2, test3,test4, test5, test6, test7, test8};

% testDataAll = {test2};
% testDataAll = {test1, test2, test3};

disp('Saving...')
save('IMUDataAll.mat','IMUDataAll')
disp('Save done.')









