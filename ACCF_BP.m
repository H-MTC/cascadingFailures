function [net,err] = ACCF_BP(ACCFM)

% samples = load('DATA_ACCFMT.mat');
% indices = max(abs(samples.ACCFM(79:end, 1:5000)),[],1) < 1000;
%   AC
% [inputs, settingsI] = mapminmax(ACCFM(1:78, :));
% [outputs, settingsO] = mapminmax(ACCFM(79:end, :));

%   DC
[inputs, settingsI] = mapminmax(ACCFM(1:58, :));
[outputs, settingsO] = mapminmax(ACCFM(59:end, :));
tr.reverse = settingsO;

% 
% inputs = inputs(:, indices);
% outputs = outputs(:,indices);

epochs = 6000;
net = bpNetConfig(200, 'trainscg', epochs, 20);
net = train(net, inputs, outputs);
% net = train(net, inputs, outputs, 'useGPU', 'yes');
tr.net = net;
save('trDC.mat','tr');

predict = sim(net, inputs);
predicted = mapminmax.reverse(predict,settingsO);
err = sum(sum(abs(predicted - outputs)));
% QQmail2me('TrainingFinished',['Training finished..., total error = ',num2str(err)]);
end

function net = bpNetConfig(numUnits, method, epochs, maxfail)

% net.layers{1}.transferFcn = 'hardlim';
% net.inputs{1}.size = size(dataTrainedX, 1);
% net.inputs{1}.processFcns = {'mapminmax','mapminmax'};
% net.outputs{2}.processFcns = {'mapminmax','mapminmax'};
net = feedforwardnet(numUnits);
net.trainFcn = method;
net.trainParam.epochs = epochs;
net.trainParam.goal	= 0;
net.trainParam.max_fail	= maxfail; %Maximum validation failures
net.trainParam.show	= 50; %Epochs between displays (NaN for no displays)
net.trainParam.showCommandLine = true; %Generate command-line output
net.trainParam.showWindow = true;
net.trainParam.time	= inf;
net.trainParam.min_grad = 1e-8;
end

