function lstmDC = dccf_lstm(dataTrain)


mu = mean(dataTrain, 2);
sig = std(dataTrain, 0, 2);

lstmDC.mu = mu;
lstmDC.sig = sig;
dataTrainStandardized = (dataTrain - repmat(mu,[1, size(dataTrain, 2)])) ./ repmat(sig,[1, size(dataTrain, 2)]);

XTrain = dataTrainStandardized(1:58, :);
YTrain = dataTrainStandardized(59:end, :);


numFeatures = 58;
numResponses = 38;
numHiddenUnits = 200;

layers = [ ...
    sequenceInputLayer(numFeatures)
    fullyConnectedLayer(numResponses)
    tanhlayer
    fullyConnectedLayer(numResponses)
    regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',5000, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',1, ...
    'ExecutionEnvironment', 'cpu');

net = trainNetwork(XTrain,YTrain,layers,options);
lstmDC.net = net;

save('lstmDC.mat','lstmDC');
end

