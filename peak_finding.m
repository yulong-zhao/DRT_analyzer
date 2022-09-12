%close all;

global pks locs


[pks_i,locs_i,w_i,p_i] = findpeaks(y);

[pks,locs,w,p] = findpeaks(y,x);

%plot(x, y); hold on;

%plot(locs, pks, 'o');


%% Plot peaks things:
% for i=1:length(pks)
%     plot([locs(i)-w(i)*0.5,locs(i)+w(i)*0.5],i*[0.1,0.1],':','Linewidth',4);
%     
%     
% end


%%



n_peaks = length(pks);

i=4;

%plot([locs(i)-w(i), locs(i)+w(i)], [0.2,0.2], '-x');


% an = annotation('doublearrow');
% [xaf,yaf] = ds2nfu(locs(i)-w(i),0.1);
% [xaf2,~] = ds2nfu(locs(i)+w(i),0.1);


%an.Position = [xaf, yaf, xaf2-xaf, 0];    
x_sep = sep_peaks(y,2);
% for i=1:size(x_sep,1)
%    figure;
%    plot(y(x_sep(i,1):x_sep(i,2) ));
%     
% end

i=1;
x1 = x(x_sep(i,1):x_sep(i,2)); 
y1 = y(x_sep(i,1):x_sep(i,2)); 

i=i+1;
x2 = x(x_sep(i,1):x_sep(i,2)); 
y2 = y(x_sep(i,1):x_sep(i,2)); 

i=i+1;
x3 = x(x_sep(i,1):x_sep(i,2)); 
y3 = y(x_sep(i,1):x_sep(i,2)); 

i=i+1;
x4 = x(x_sep(i,1):x_sep(i,2)); 
y4 = y(x_sep(i,1):x_sep(i,2)); 

i=i+1;
x5 = x(x_sep(i,1):x_sep(i,2)); 
y5 = y(x_sep(i,1):x_sep(i,2)); 

i=i+1;
x6 = x(x_sep(i,1):x_sep(i,2)); 
y6 = y(x_sep(i,1):x_sep(i,2)); 

% %%
% for i=1:length(locs_i)
%    figure;
%    plot(y(locs_i(i)-floor(w_i(i)):locs_i(i)+floor(w_i(i))));
%     
% end

if(cand_info.mode == "single")
    x_try = x(x_sep(i_try,1):x_sep(i_try,2));
    y_try = y(x_sep(i_try,1):x_sep(i_try,2));
    
elseif(cand_info.mode == "first")
    x_try = x_c(1:369);
    y_try = y_c(1:369);
else
    x_try = x_c;
    y_try = y_c;
end

peaks_sel = x_sep(i_try,1)<= locs_i &  locs_i<= x_sep(i_try,2);
locs_try = locs ( peaks_sel);
w_try = w(peaks_sel);
n_locs_try = length(locs_try);