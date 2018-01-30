function clean_up_workspace()
    clear;
    clc;
    delete(findall(0,'Type','figure'));
end
