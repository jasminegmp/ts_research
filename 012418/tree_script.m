load('merge_info_end.mat')
load('merge_info_start.mat')
load('merge_info_count.mat')

[count,count_idx] = sort(merge_info.merge_count,'ascend');

for i = 1:length(merge_info.start)
    
end