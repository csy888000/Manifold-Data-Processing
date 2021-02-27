clear
close all

addpath gpdm gpdm/netlab gpdm/gpdm gpdm/gplvm gpdm/util

format long

global USE_GAMMA_PRIOR  % gamma prior for dynamics, only works with RBF kernel
global GAMMA_ALPHA % defines shape of the gamma prior
global USE_LAWRENCE % fix dynamics HPs, as Lawrence suggested (use with thetad = [0.2 0.01 1e6];) 
global FIX_HP % fix all HPs
global MARGINAL_W % marginalize over W while learning X
global MARGINAL_DW % marginalize over scale in dynamics while learning X
global LEARN_SCALE % use different scales for different output dimensions
global REMOVE_REDUNDANT_SCALE % let W absorb the overall scale of reconstruction
global W_VARIANCE % kappa^2 in the paper, not really the variance though
global M_CONST % M value in Jack's master's thesis
global BALANCE % Constant in front of dynamics term, set to D/q for the B-GPDM
global SUBSET_SIZE % Number of data to select for EM, set -1 for all data. 
global USE_OLD_MISSING_DATA

M_CONST = 1; 
REMOVE_REDUNDANT_SCALE = 1;
LEARN_SCALE = 1; 
MARGINAL_W = 0; 
MARGINAL_DW = 0; 
W_VARIANCE = 1e6; 
FIX_HP = 0; 
USE_GAMMA_PRIOR = 0; 
GAMMA_ALPHA = [5 10 2.5]; 
USE_LAWRENCE = 0;
BALANCE = 1;
SUBSET_SIZE = -1; 


opt = foptions;
opt(1) = 1;
opt(9) = 0;
if MARGINAL_W == 1
    opt(14) = 100; % total number of iterations
    extItr = 1; 
else
    opt(14) = 10; % rescaling every 10 iterations
    extItr = 100; % do extItr*opt(14) iterations in total
end  

% modelType(1) : input of dynamics
%   0 => [x_t, x_{t-1}]
%   1 => [x_t, x_t - x_{t-1}]
%   2 => [x_t]
% modelType(2) : output of dynamics 
%   0 => x_{t+1} 
%   1 => x_{t+1} - x_t
% modelType(3) : kernel type
%   0 => RBF kernel with weighted dimensions, use with input 0 or 1
%   1 => RBF kernel 
%   2 => Linear kernel
%   3 => weighted Linear kernel + RBF kernel with weighted dimensions, use with
%   input 0 or 1
%   4 => weighted linear kernel
%   5 => linear + RBF

% Learn single walker model from lin+rbf kernel.
modelType = [2 0 5]; 


%%







%%
load MatSave/walkAll.mat

num = 8;
labelNum = 4;
% angleLenAll = [115 131 130 134];
% angleLenAll = [70 40 30 256];
angleLenAll = [120 124 128 128];



heelLeftValueSave = walkAll{num}.heelLeftValueSave;
heelRightValueSave = walkAll{num}.heelRightValueSave;
walkSeg = walkAll{num}.walkSeg;
imuSeg = walkAll{num}.imuSeg;
jointTh = walkAll{num}.jointTh;
imuTrial = walkAll{num}.imuTrial;
label_segment = walkAll{num}.label_segment;



k = 1;
for i = labelNum:4:length(walkSeg)
    heelLeftV = heelLeftValueSave{i};
%     for j = 1:(length(heelLeftV)-1)
    if length(heelLeftV) > 2
        j = 1;    
        heelLeftStep{k} = heelLeftV(j):heelLeftV(j+1);
        jointAngle{k} = jointTh(heelLeftStep{k},:);
        imuData{k} = imuTrial(heelLeftStep{k},:);
        k = k+1;
    end
end



figure(21)
clf
% colorX = {'b','k','r','g','m','c'};
for i = 1:length(jointAngle)
    disp(['Current is NUMBER: ', num2str(i)])
    for j = 1:size(jointAngle{1},2)
        subplot(3,4,j)
        hold on
        plot(jointAngle{i}(:,j))
    end
end

%%
lenA = 0;
for i = 1:length(jointAngle)
    lenA = lenA + length(jointAngle{i})/length(jointAngle);
end
% angleLen = round(lenA);

angleLen = angleLenAll(labelNum);
% angleLen = 115;
jointAngleMean = zeros(angleLen,size(jointAngle{1},2));
for i = 1:length(jointAngle)
    for j = 1:size(jointAngle{1},2)
        angleOne = jointAngle{i}(:,j);   
        angleOne = interp1(1:length(angleOne), angleOne,...
            0:length(angleOne)/(angleLen-1):length(angleOne),'spline')';
        jointAngle_interp{i}(:,j) = angleOne;
        jointAngleMean(:,j) = jointAngleMean(:,j) + angleOne./length(jointAngle);
    end
end

phase = 0:(1/(angleLen-1)):1;
jointAngleMeanAcross = [phase' jointAngleMean];
writematrix(jointAngleMeanAcross, 'jointAngleMeanAcross.txt')

%%
figure(22)
clf
hold on
for j = 1:size(jointAngle{1},2)
    subplot(3,4,j)
    hold on
    plot(jointAngleMean(:,j))
end

jointAngleOne{1} = jointAngleMean;


%%
walkX = cell(1,length(jointAngleOne));
walkX_pred = cell(1,length(jointAngleOne));
for j = 1:length(jointAngleOne)
    Y = jointAngleOne{j};
    segments = 1;

    missing = [];
    N = size(Y, 1); D = size(Y, 2);
    q = 3; % dimensionality of latent space

    % PCA
    X = zeros(N, q);
    refY = Y; meanData = mean(Y);
    Y = Y - repmat(meanData, N, 1);
    [v, u] = pcaGP(Y);
    [coeff,score,latent] = pca(Y);
    v(v<0)=0;
    X = Y*u(:, 1:q)*diag(1./sqrt(v(1:q)));

    % initialize hyperparameters
    theta = [1 1 exp(1)];
    thetad = [0.9 1 0.1 exp(1)];
    w = ones(D,1);

    [X,theta,thetad,w] = gpdmfitFull(X, Y, w, segments, theta, thetad, opt, ... 
         extItr, modelType, missing);

    % load example_model
    [K,invK] = computeKernel(X, theta);
    [Xin,Xout] = priorIO(X, segments, modelType);
    [Kd,invKd] = computePriorKernel(Xin, thetad, modelType(3));
    simSteps = 128;
%     simSteps = length(jointAngleOne{1});
    % starts at ened of training sequence;
    simStart = [X(segments(1)+1,:) X(end,:)]; %  inputs 2 points in case using 2nd order model
    [X_pred,XRand_pred] = simulatedynamics(X, segments, thetad, invKd, simSteps, simStart, modelType);


    
    walku{j} = u;
    walkv{j} = v;
    walkY{j} = Y;
    walkX{j} = X;
    walkX_pred{j} = X_pred;
    
end





%%


% segments = 1;

figure(10)
clf;
hold on;

colorX = {'b','k','r','g','m','c'};
colorXpred = {'b','k','r','g','m','c'};
shapeU = {'d','s','x','+','*','^'};

walk_pred_rot = walkX_pred;
for k = 1%[2 4 5 6]%1:length(jointAngleOne)
    
    X = walkX{k};
    X_pred = walkX_pred{k};

%     X(:,[1 2 3]) = X(:,[2 1 3]);
%     X(:,3) = -X(:,3);
%     
%     X_pred(:,[1 2 3]) = X_pred(:,[2 1 3]);
%     X_pred(:,3) = -X_pred(:,3);

    walk_pred_rot{k} = X_pred;
    plotseries(X_pred, [1], colorXpred{k});
    
    hold on

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    view(9,8)

end
plot3(X_pred(1,1),X_pred(1,2),X_pred(1,3),'ro')
plot3(X_pred(end,1),X_pred(end,2),X_pred(end,3),'gd')



figure(23)
plot3(X_pred(:,1),X_pred(:,2),X_pred(:,3))
axis equal

manifoldData.Xpred = X_pred;
manifoldData.jointAngleMean = jointAngleMean;

save('manifoldData.mat','manifoldData')

%%
k = 1;
for i = labelNum:4:length(walkSeg)
    heelLeftV = heelLeftValueSave{i};
    for j = 1:(length(heelLeftV)-1)
        heelLeftStepAll{k} = heelLeftV(j):heelLeftV(j+1);
        jointAngleAll{k} = jointTh(heelLeftStepAll{k},:);
        imuDataAll{k} = imuTrial(heelLeftStepAll{k},:);
        k = k+1;
    end
end
    

for i = 1:length(jointAngleAll)
    lenData = size(jointAngleAll{i},1);
    for j = 1:3
        predData = X_pred(:,j);
        X_pred_interp = interp1(1:length(predData), predData,...
            0:length(predData)/(lenData-1):length(predData),'spline')';
        X_pred_all{i}(:,j) = X_pred_interp;
    end
    X_propagation{i} = 0:(1/(lenData-1)):1;
end


% figure(24)
% hold on
% for i = 1:length(X_pred_all)
%     plot3(X_pred_all{i}(:,1),X_pred_all{i}(:,2),X_pred_all{i}(:,3))
% end

trainIMU = [];
trainXpred = [];
trainPropagation = [];
trainAngle = [];
for i = 1:length(X_pred_all)
    trainIMU = [trainIMU;imuDataAll{i}];
    trainXpred = [trainXpred;X_pred_all{i}];
    trainPropagation = [trainPropagation;X_propagation{i}'];
    trainAngle = [trainAngle;jointAngleAll{i}];
end


labelData = ones(1,length(trainIMU))';
trainData_walkAcross = [trainIMU trainXpred trainPropagation trainAngle labelData];

save('trainData_walkAcross.mat','trainData_walkAcross')
writematrix(trainData_walkAcross,'trainData_walkAcross.csv') 




%%

function rotMat3 = rotAngle(x, y, z)

rotMatX = [1 0 0; 0 cosd(x) sind(x); 0 -sind(x) cosd(x)];
rotMatY = [cosd(y) 0 -sind(y); 0 1 0; sind(y) 0 cosd(y)];
rotMatZ = [cosd(z) sind(z) 0; -sind(z) cosd(z) 0; 0 0 1];

rotMat3 = rotMatX*rotMatY*rotMatZ;
end


