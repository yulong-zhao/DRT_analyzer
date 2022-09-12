dx_max = diff(data(end-1:end,1));
dx_min = dx_max/15; 


x_data = data(1:end,1); % Beginning is zero. 
y_data = data(1:end,2);

% bnd_zero = find(diff(y_data<1e-9));
% bnd_zero2 = bnd_zero;
% bnd_zero2(1) = bnd_zero2(1) - 5;
% bnd_zero2(end) = bnd_zero2(2) + 5;
% 
% 
% x_data = x_data(bnd_zero2(1):bnd_zero2(end));
% y_data = y_data(bnd_zero2(1):bnd_zero2(end));

scaling.x.min = min(x_data);
x_data_n = (x_data - min(x_data));

scaling.x.max = max(x_data_n);
x_data_n = x_data_n/max(x_data_n);

scaling.y.min = min(y_data);
y_data_n = (y_data - min(y_data));

scaling.y.max = max(y_data_n);
y_data_n = y_data_n/max(y_data_n);


r_data_n = vecnorm(diff([x_data_n,y_data_n]),2,2);


%(logspace(0,log10(x_data(end)-x_data(1)+1),n_data)-1)+x_data(1); %linspace(data(1,1), data(end,1),1500);%    
%(); linspace(data(1,1), data(end,1),1500); % (logspace(0,log10(data(end,1)-data(1,1)+1),n_data)-1)+data(1,1); %linspace(x_data(1), x_data(end), n_data); 
n_data = 3000;
x =  linspace(x_data(1), x_data(end), n_data); 
x = x(:);
%vec(linspace(data(1,1), data(end,1),1000));% data(:,1);
y_orj = interp1(x_data, y_data, x,'makima');% data(:,2);

y = y_orj /max(y_orj);

% %%
% figure;
% plot(data(:,1),data(:,2)); hold on;
% plot(x,y,'-.');

% x_c = x(430:end);
% y_c = y(430:end);

x_try = x;
y_try = y;
