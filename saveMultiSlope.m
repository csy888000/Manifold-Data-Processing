clear
close all

%%
manifoldData{1} = load('manifoldData_5deg.mat');
manifoldData{2} = load('manifoldData_10deg.mat');
manifoldData{3} = load('manifoldData_15deg.mat');
manifoldData{4} = load('manifoldData_20deg.mat');




%%
figure(1)
clf
hold on;
colorX = {'b','k','r','m','c','g'};
colorXpred = {'b','k','r','m','c','g'};
shapeU = {'d','s','x','+','*','^'};

for i = 1:4
    X_pred = manifoldData{i}.manifoldData.Xpred;

    
    
    if i == 1
        X_pred = X_pred*rotAngle(0,0,0);
    end
    
    if i == 2
        X_pred = X_pred*rotAngle(-35,0,0);
    end

    
    if i == 3
        X_pred = X_pred*rotAngle(-40,20,0);
    end
    
    if i == 4
        X_pred(:,2) = -X_pred(:,2);
        X_pred = X_pred*rotAngle(-50,0,0);
    end
    
    X_predNew{i} = X_pred;
    plotseries(X_pred, [1], colorXpred{i});
    hold on
    

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    view(49,18)
    
    plot3(X_pred(1,1),X_pred(1,2),X_pred(1,3),'ro')
    plot3(X_pred(end,1),X_pred(end,2),X_pred(end,3),'gd')
end

%%
figure(2)
clf
for i = 1:4
    X_pred = X_predNew{i};
    for k = 1:3
        subplot(3,1,k)
        hold on
        plot(X_pred(:,k), colorXpred{i})
    end
end
        
    
%%
figure(3)
clf
for i = 1:4
    jointAngleMean = manifoldData{i}.manifoldData.jointAngleMean;
    for k = 1:12
        subplot(3,4,k)
        hold on
        plot(jointAngleMean(:,k), colorXpred{i})
    end
end








%%


function rotMat3 = rotAngle(x, y, z)

rotMatX = [1 0 0; 0 cosd(x) sind(x); 0 -sind(x) cosd(x)];
rotMatY = [cosd(y) 0 -sind(y); 0 1 0; sind(y) 0 cosd(y)];
rotMatZ = [cosd(z) sind(z) 0; -sind(z) cosd(z) 0; 0 0 1];

rotMat3 = rotMatX*rotMatY*rotMatZ;
end