function anchor = getAnchor(fea,numofsample)

[n, ~] = size(fea);
data_cell = {fea};
%% splitting
while 1
    ball_number_old = length(data_cell);
    data_cell = division(data_cell);
    ball_number_new = length(data_cell);
    if ball_number_new == ball_number_old 
       break
    end
end
m=length(data_cell);
last_center={};
last_radius={};
temp_center={};
%% filter
MaxAnchorNum=ceil(numofsample/10);
currentLength = 1;
while m>MaxAnchorNum
    for h=1:m 
        if size(data_cell{h},1)==currentLength
            continue;
        end
        temp_center{end+1,1}=data_cell{h};
    end   
    if length(temp_center) == 0
        break;
    else
        m=length(temp_center);
    end
    data_cell = temp_center;
    temp_center={};
    currentLength=currentLength+1;
end

for h=1:m
    data = cell2mat(data_cell{h});
    if size(data,1)==1   
%         continue; 
        last_center{end+1,1}=data;
    else
        last_center{end+1,1}=mean(data);
        last_radius{end+1,1}=get_radius(data);
    end
end
anchor=cell2mat(last_center);
last_radius=cell2mat(last_radius);
end

%% 分裂球
function gb_newcell = division(hb_cell)
gb_newcell={};
i=1;
 for j =1:length(hb_cell)
     [m,~]=size(hb_cell{j});
     if m>8
            [ball_1, ball_2] = split_ball(hb_cell{j});
            %split by WDV
            parent_dm = get_density_volume(hb_cell{j});
            child_1_dm = get_density_volume(ball_1);
            child_2_dm = get_density_volume(ball_2);
            w=size(ball_1,1)+size(ball_2,1);
            w1=size(ball_1,1)/w;
            w2=size(ball_2,1)/w;
            w_child_dm = (w1 * child_1_dm + w2 * child_2_dm);  % WDM
            t1 = (child_1_dm > parent_dm) & (child_2_dm > parent_dm);
            t2 = (w_child_dm > parent_dm);  
            t3 = (size(ball_1, 1) > 0) && (size(ball_2, 1) > 0);
            if t2
                gb_newcell{i,1} = ball_1;
                gb_newcell{i+1,1} = ball_2;
                i = i+ 2;
            else
                gb_newcell{i,1} = hb_cell{j};
                i = i+1;
            end  
     else
            gb_newcell{i,1} = hb_cell{j};
            i = i+1;  
     end
 end
end

%% 2-split
function [ball_1, ball_2] = split_ball(data)
    if iscell(data)
        data = cell2mat(data);
    end

    % Step 1: compute center
    c = mean(data);

    % Step 2: p is farest from c
    distances_to_c = pdist2(data, c);
    [~, idx_p] = max(distances_to_c);
    p = data(idx_p, :);

    % Step 3: q is farest from p
    distances_to_p = pdist2(data, p);
    [~, idx_q] = max(distances_to_p);
    q = data(idx_q, :);

    % Step 4: compute new center c1 和 c2
    c1 = (c + p) / 2;
    c2 = (c + q) / 2;

    % Step 5: assign data to subBall according to distance
    ball_1 = {};
    ball_2 = {};
    i = 1;
    k = 1;
    
    for j = 1:size(data, 1)
        dist_to_c1 = norm(data(j, :) - c1);
        dist_to_c2 = norm(data(j, :) - c2);
        
        if dist_to_c1 < dist_to_c2
            ball_1{i, 1} = data(j, :);
            i = i + 1;
        else
            ball_2{k, 1} = data(j, :);
            k = k + 1;
        end
    end
end

%% getRadius

function radius=get_radius(data)
[num,~]=size(data);
center=mean(data);
diffMat=repmat(center,num,1);
sqDiffMat=(diffMat-data).^2;
sqDistances =sum(sqDiffMat);
distances=sqrt(sqDistances);
radius=max(distances);
end
%% get DV
function density_volume = get_density_volume(gb)
    if iscell(gb)
        gb = cell2mat(gb);
    end
   
    num = size(gb, 1);
    density_volume = 0;
    
    if num > 0
        center = mean(gb, 1);
        diffMat = repmat(center, num, 1) - gb;
        sqDiffMat = diffMat .^ 2;
        sqDistances = sum(sqDiffMat, 2);
        distances = sqrt(sqDistances);
        sum_radius = sum(distances);
        
        if sum_radius ~= 0
            density_volume = num / sum_radius;
        else
            density_volume = num; % 当半径和为0时，直接用num表示密度体积
        end
    end
end

