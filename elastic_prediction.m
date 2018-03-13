function out = elastic_prediction(f, t, model, option)
% ELASTIC_PREDICTION Elastic Functional Regression Prediction
% -------------------------------------------------------------------------
% This function identifies a logistic regression model with
% phase-variablity using elastic methods
%
% Usage:  out = elastic_prediction(f, t, model)
%         out = elastic_prediction(f, t, model, option)
%
% Input:
% f (M,N): matrix defining N functions of M samples
% t : time vector of length M
% model: model structure from elastic_regression, elastic_logistic, or
% elastic_mlogistic
%
% default options
% option.y = []; % optional truth
% option.smooth = 0; % smooth data using standard box filter
% option.sparam = 25; % number of times to run filter
%
% Output:
% structure with fields:
% y_pred: predicted value or probability (depends on model type)
% SSE: if linear regression
% y_labels: assigned labels (logistic and mlogistic)
% PC: probability of classification (logistic and mlogistic)

if nargin < 4
    option.y = [];
    option.smooth = 0;
    option.sparam = 25;
end

if option.smooth == 1
    f = smooth_data(f, option.sparam);
end

q = f_to_srvf(f,t);
M = length(t);
n = size(q,2);
SSE = NaN;
PC = NaN;
y_labels = [];

if strcmpi(model.type, 'linear') || strcmpi(model.type,'logistic')
    y_pred = zeros(n,1);
elseif strcmpi(model.type, 'mlogistic')
    m = model.n_classes;
    y_pred = zeros(n,m);
else
    error('unkown model type')
end

for ii = 1:n
    difference = model.q - repmat(q(:,ii),1,size(model.q,2));
    dist = sum(abs(difference).^2).^(1/2);
    [~, argmin] = min(dist);
    gamdev = gradient(model.gamma(:,argmin), 1/(M-1));
    q_tmp = interp1(t, q(:,ii), (t(end)-t(1)).*model.gamma(:,argmin) ...
        + t(1))'.*sqrt(gamdev');
    
    switch model.type
        case 'linear'
            y_pred(ii) = model.alpha + trapz(t, q_tmp' .* model.beta);
        case 'logistic'
            y_pred(ii) = model.alpha + trapz(t, q_tmp' .* model.beta);
        case 'mlogistic'
            for jj = 1:m
                y_pred(ii,jj) = model.alpha(jj) + trapz(t, q_tmp' .* model.beta(:, jj));
            end
    end
end

if isempty(option.y)
    switch model.type
        case 'linear'
            SSE = NaN;
        case 'logistic'
            y_pred = phi(y_pred);
            y_labels = ones(1,n);
            y_labels(y_pred < 0.5) = -1;
            PC = NaN;
        case 'mlogistic'
            y_pred = phi(reshape(y_pred,1,n*m));
            y_pred = reshape(y_pred,n,m);
            [~, y_labels] = max(y_pred,[],2);
            PC = NaN;
    end
else
    switch model.type
        case 'linear'
            SSE = sum((option.y-y_pred).^2);
        case 'logistic'
            y_pred = phi(y_pred);
            y_labels = ones(1,n);
            y_labels(y_pred < 0.5) = -1;
            TP = sum(option.y(y_labels == 1) == 1);
            FP = sum(option.y(y_labels == -1) == 1);
            TN = sum(option.y(y_labels == -1) == -1);
            FN = sum(option.y(y_labels == 1) == -1);
            PC = (TP+TN)/(TP+FP+FN+TN);
        case 'mlogistic'
            y_pred = phi(reshape(y_pred,1,n*m));
            y_pred = reshape(y_pred,n,m);
            [~, y_labels] = max(y_pred,[],2);
            PC = zeros(1,m);
            cls_set = 1:m;
            for ii = 1:m
                cls_sub = setdiff(cls_set,ii);
                TP = sum(option.y(y_labels == ii) == ii);
                FP = sum(option.y(ismember(y_labels,cls_sub)) == ii);
                TN = sum(option.y(ismember(y_labels,cls_sub)) == ...
                    y_labels(ismember(y_labels,cls_sub)));
                FN = sum(ismember(option.y(y_labels==ii), cls_sub));
                PC(ii) = (TP+TN)/(TP+FP+FN+TN);
            end
            PC = sum(option.y == y_labels)./length(y_labels);
    end
end
out.y_pred = y_pred;
out.y_labels = y_labels;
out.SSE = SSE;
out.PC = PC;
end

%% Helper functions

function out = phi(t)
% calculates logisitc function, returns 1/(1+exp(-t))
idx = t > 0;
out = zeros(size(t));
out(idx) = 1./(1+exp(-t(idx)));
exp_t = exp(t(~idx));
out(~idx) = exp_t ./ (1+exp_t);
end
