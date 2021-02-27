clc
clear



subjects{1}.age = 26;
subjects{1}.height = 165;
subjects{1}.mass = 64;

subjects{2}.age = 32;
subjects{2}.height = 181;
subjects{2}.mass = 76;

subjects{3}.age = 31;
subjects{3}.height = 172;
subjects{3}.mass = 74;


subjects{4}.age = 32;
subjects{4}.height = 170;
subjects{4}.mass = 79;


age_array = zeros(1,length(subjects));
height_array = zeros(1,length(subjects));
mass_array = zeros(1,length(subjects));
for i = 1:length(subjects)
    age_array(i) = subjects{i}.age;
    height_array(i) = subjects{i}.height;
    mass_array(i) = subjects{i}.mass;
end

age_average = mean(age_array);
age_variance = std(age_array);
height_average = mean(height_array);
height_variance = std(height_array);
mass_average = mean(mass_array);
mass_variance = std(mass_array);


disp(['age_average is ', num2str(age_average)])
disp(['age_variance is ', num2str(age_variance)])
disp(['height_average is ', num2str(height_average)])
disp(['height_variance is ', num2str(height_variance)])
disp(['mass_average is ', num2str(mass_average)])
disp(['mass_variance is ', num2str(mass_variance)])

