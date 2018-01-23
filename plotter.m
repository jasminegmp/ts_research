load('tree.mat');

for level = 1:21
    figure;
    hold on;
    title(level)
    for idx = 1:length(tree)
        if (tree(idx).level == level) || (tree(idx).level == level+1)
            time = tree(idx).index(1,1):1:tree(idx).index(1,2);
            if (tree(idx).parent_idx == 0) && (tree(idx).level == level+1)
                plot(time, tree(idx).data, 'Color', [1 0 0], 'LineWidth', 1);
            elseif (tree(idx).level == level) && (tree(idx).parent_idx ~= 0)
                plot(time, tree(idx).data, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7);
            end
        end
    end
end