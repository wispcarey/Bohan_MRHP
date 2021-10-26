clear;
clc;

%% 4/18/2021

A = readtable('covid19cases_test.csv');

B = readtable('uscities.csv');
% Btext_info = B.textdata;
% Bvalue_info = B.data;

loc_name = table2cell(unique(A(:,2)));

hasMatch = cellfun(@(x) strcmp(x,'CA'), table2cell(B(:,'state_id')));
loc = B(hasMatch,{'county_name','lat','lng'});

loc_info = zeros(size(loc_name,1),2);

for i = 1:size(loc_name,1)
    ind = cellfun(@(x) strcmp(x,loc_name(i)), table2cell(loc(:,'county_name')));
    if any(ind)
        loc_info(i,1) = mean(cellfun(@str2num, table2cell(loc(ind, 'lat'))));
        loc_info(i,2) = mean(cellfun(@str2num, table2cell(loc(ind, 'lng'))));
    end
end

% loc_name = loc_name(loc_info(:,1)~=0);
% loc_info = loc_info(loc_info(:,1)~=0,:);

dist_to_la = sum((loc_info - loc_info(cellfun(@(x) strcmp(x, 'Los Angeles'), loc_name), :)).^2,2);
[~,I] = sort(dist_to_la);

%% choose K closest locations to LA 
K = 10;
top_loc_name = loc_name(I(1:K));
selected_inds = cellfun(@(x) any(cellfun(@(y) strcmp(y, x), top_loc_name)), table2cell(A(:,'area')));
selected_A = A(selected_inds, :);

%% timeline
T = cellfun(@(x) daysact(string(table2cell(selected_A(end,'date'))),x), table2cell(selected_A(:,'date')));
save timeline T

%% location to index

Area_ind = cellfun(@(x) replace_cell_ind(x, top_loc_name), table2cell(selected_A(:,'area')));
save Area_ind Area_ind

%% use reported_cases
D_start = 1;
D_end = 435;
N_size = 1;
T1 = flip(T(T>=D_start& T<=D_end))-D_start;
Area_ind = flip(Area_ind(T>=D_start& T<=D_end));
RC = cellfun(@(x)x, table2cell(selected_A(:,'reported_cases')));
RC = flip(RC(T>=D_start& T<=D_end));
RC = ceil(RC/N_size);
RC(RC<0) = 0;

timeline = zeros(sum(RC),1);
types = zeros(sum(RC),1);

k = 1;
for i = 1:length(RC)
    if RC(i)~=0
        timeline(k:k+RC(i)-1) = T1(i);
        types(k:k+RC(i)-1) = Area_ind(i);
        k = k + RC(i);
    end
end

%% cumulative cases
Cum_cases = zeros(K,D_end-D_start+2);
for i = 1:sum(RC)
    Cum_cases(types(i), timeline(i)+2:end) = Cum_cases(types(i), timeline(i)+2:end) + 1;
end
    
save original_data timeline types Cum_cases


function I = replace_cell_ind(x, C)
    I = 0;
    for i = 1:length(C)
        if strcmp(x, C(i))
            I = i;
            break;
        end
    end
end