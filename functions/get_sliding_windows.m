function groups = get_sliding_windows(N,window_length,step_length)

window_num = floor((N - window_length) / step_length) + 1;

for i = 1:window_num
    groups{i,1} = [1:window_length] + (i-1)*step_length;
end

