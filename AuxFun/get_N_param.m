function n_param = get_N_param(function_handle)

text = func2str(function_handle); 
i = find('1'<=text & text <='9',1,'last');

n_param = str2num( text(i));

end