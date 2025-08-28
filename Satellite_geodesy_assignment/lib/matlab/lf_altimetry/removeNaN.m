
function [d a] = removeNaN(data)

%% PURPOSE:
% remove NaN values from a vector

%% INPUT:
% data: vector with NaN

%% OUTPUT:
% d: vector without NaN elements
% a: number of removed values


%%
index = find(~isnan(data)); 
d=data(index);
a=size(data,1)-size(d,1);
