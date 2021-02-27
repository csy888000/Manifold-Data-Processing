function [Vicon_filled, ViconObj] = fillViconGapsFullbody(dataMarkers)

dataMarkers(dataMarkers == 0) = NaN;

% markerNames = {'LASI', 'RASI', 'LPSI', 'RPSI', 'LTHI', 'LKNE', 'LTIB', 'LANK',...
%     'LHEE', 'LTOE', 'RTHI', 'RKNE', 'RTIB', 'RANK', 'RHEE', 'RTOE'};


markerNames = {'LFHD','RFHD','LBHD','RBHD',...
    'C7','T10','CLAV','STRN','RBAK',...
    'LSHO','LUPA','LELB','LFRM','LWRA','LWRB','LFIN',...
    'RSHO','RUPA','RELB','RFRM','RWRA','RWRB','RFIN',...
    'LASI','RASI','LPSI','RPSI',...
    'LTHI','LKNE','LTIB','LANK','LHEE','LTOE',...
    'RTHI','RKNE','RTIB','RANK','RHEE','RTOE'};



for i = 1:length(markerNames)
    oneMarkerData = dataMarkers(:,i*3-2:i*3);
    eval([(markerNames{i}), ' = oneMarkerData;'])
    eval([markerNames{i} '= fillgaps(', markerNames{i}, ');']);
%     eval(['LFHD = oneMarkerData'])
    eval(['ViconObj.',markerNames{i}, '=', markerNames{i}, ';']);
end


Vicon_filled = [LFHD RFHD LBHD RBHD C7 T10 CLAV STRN RBAK,...
    LSHO LUPA LELB LFRM LWRA LWRB LFIN RSHO RUPA RELB RFRM RWRA RWRB RFIN,...
    LASI RASI LPSI RPSI LTHI LKNE LTIB LANK LHEE LTOE RTHI RKNE RTIB RANK RHEE RTOE];


% LASI = dataMarkers(:,1:3);
% RASI = dataMarkers(:,4:6);
% LPSI = dataMarkers(:,7:9);
% RPSI = dataMarkers(:,10:12);
% 
% LTHI = dataMarkers(:,13:15);
% LKNE = dataMarkers(:,16:18);
% LTIB = dataMarkers(:,19:21);
% LANK = dataMarkers(:,22:24);
% LHEE = dataMarkers(:,25:27);
% LTOE = dataMarkers(:,28:30);
% 
% RTHI = dataMarkers(:,31:33);
% RKNE = dataMarkers(:,34:36);
% RTIB = dataMarkers(:,37:39);
% RANK = dataMarkers(:,40:42);
% RHEE = dataMarkers(:,43:45);
% RTOE = dataMarkers(:,46:48);



%%
% LASI = fillgaps(LASI);
% RASI = fillgaps(RASI);
% LPSI = fillgaps(LPSI);
% RPSI = fillgaps(RPSI);
% 
% LTHI = fillgaps(LTHI);
% LKNE = fillgaps(LKNE);
% LTIB = fillgaps(LTIB);
% LANK = fillgaps(LANK);
% LHEE = fillgaps(LHEE);
% LTOE = fillgaps(LTOE);
% 
% RTHI = fillgaps(RTHI);
% RKNE = fillgaps(RKNE);
% RTIB = fillgaps(RTIB);
% RANK = fillgaps(RANK);
% RHEE = fillgaps(RHEE);
% RTOE = fillgaps(RTOE);


% Vicon_filled = [LASI RASI LPSI RPSI LTHI LKNE LTIB LANK LHEE LTOE RTHI RKNE RTIB RANK RHEE RTOE];

% ViconObj.LASI = LASI;
% ViconObj.RASI = RASI;
% ViconObj.LPSI = LPSI;
% ViconObj.RPSI = RPSI;
% 
% ViconObj.LTHI = LTHI;
% ViconObj.LKNE = LKNE;
% ViconObj.LTIB = LTIB;
% ViconObj.LANK = LANK;
% ViconObj.LHEE = LHEE;
% ViconObj.LTOE = LTOE;
% 
% ViconObj.RTHI = RTHI;
% ViconObj.RKNE = RKNE;
% ViconObj.RTIB = RTIB;
% ViconObj.RANK = RANK;
% ViconObj.RHEE = RHEE;
% ViconObj.RTOE = RTOE;



end