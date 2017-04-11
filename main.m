close all; clear;
clc;

%% Data Handling:

trainFolder ='YaleSubsets/1';
trainImages = dir(strcat(trainFolder,'/*.pgm'));
trainLabels = zeros(1,length(trainImages));
for k = 1:length(trainImages)
    tmp = trainImages(k).name;
    trainLabels(k) = str2double(tmp(6:7));
end


testFolder ='YaleSubsets/5';
testImages = dir(strcat(testFolder,'/*.pgm'));
testLabels = zeros(1,length(testImages));
for k = 1:length(testImages)
    tmp = testImages(k).name;
    testLabels(k) = str2double(tmp(6:7));
end


%% Feature Extraction:

% Image parameters
r=128; c=128;
num_pixels = r*c;

% construct Gabor filterbank:
num_scales = 5;
num_orientations = 8;
filter_bank = construct_Gabor_filters_PhD(num_orientations, num_scales, [r c]);

% Histogram parameters
sub_region_size = [16 16];
num_bins = 16;

num_subregions = (r/sub_region_size(1))*(c/sub_region_size(2));

trainFeats = zeros(length(trainImages),90*num_bins*num_subregions);
for k = 1:length(trainImages);
    img = imread(strcat(trainFolder,'/',trainImages(k).name));
    img = imresize(img,[128 128],'bilinear');
    filtered_image = filter_image_with_Gabor_bank_PhD(img,filter_bank,1);
    HGPP = compute_HGPP( filtered_image, r, c, num_pixels, num_scales, num_orientations, sub_region_size, num_bins );
    trainFeats(k,:) = HGPP';
end


testFeats = zeros(length(trainImages),90*num_bins*num_subregions);
for k = 1:length(testImages);
    img = imread(strcat(testFolder,'/',testImages(k).name));
    img = imresize(img,[128 128],'bilinear');
    filtered_image = filter_image_with_Gabor_bank_PhD(img,filter_bank,1);
    HGPP = compute_HGPP( filtered_image, r, c, num_pixels, num_scales, num_orientations, sub_region_size, num_bins );
    testFeats(k,:) = HGPP';
end


%% Classification

%results = nn_classification_PhD(modelLDA.train, trainLabels,testFeatsLDA , testLabels, size(testFeatsLDA,1), 'mahcos');
scores = zeros(length(testImages),length(trainImages));

for i=1:length(testImages)
    for j=1:length(trainImages)
        hist_train = trainFeats(j,:);
        hist_test = testFeats(i,:);
        scores(i,j) = histogram_intersection(hist_train',hist_test',num_bins);
    end
end

[~,I] = max(scores,[],2);

%output = evaluate_results_PhD(results,'ID');

labels_predict = trainLabels(I);
accuracy = sum(testLabels==labels_predict)/length(testLabels)


















