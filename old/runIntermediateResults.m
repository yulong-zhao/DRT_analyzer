clear variables; close all;

myFolder = 'intermediate_results/';
for i_try =1:6

new_code_2021_04_03;

save([myFolder, 'res2021_04_04_i_try_', num2str(i_try), '.mat'], 'x_opt', 'R_opt', 'A_opt');

end