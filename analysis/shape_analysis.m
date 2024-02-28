clear;clc;close all; 
load('Z:\Jessica\tracking_code\160523\data_1\1_6_1_mask_tracking_varied_model.mat')
load('Z:\Jessica\tracking_code\160523\label_matrix\data\1_6_1\labels\10_filtered_labels_consistent_1_6_1.mat')
common_val = nonzeros(common_val);
mask_track_x = mask_track(:,:,size(mask_track,3)-4);
mask_track_y = mask_track(:,:,size(mask_track,3)-3);
nuc_track_x = mask_track(:,:,1);
nuc_track_y = mask_track(:,:,2);
mask_area = mask_track(:,:,13);
mask_perimeter = mask_track(:,:,18);
v = mask_track(:,:,6);
coord_front = mask_track(:,:,7);
coord_front(isnan(coord_front)) = 0;
avg_front = mean(coord_front,2);

coord_lat = mask_track(:,:,8);
coord_lat(isnan(coord_lat)) = 0;
avg_lat = mean(coord_lat,2);
coord_back = mask_track(:,:,9);
coord_back(isnan(coord_back)) = 0;
avg_back = mean(coord_back,2);
ecc = mask_track(:,:,17);
circ = mask_track(:,:,18);
ecc(isnan(ecc)) = 0;
circ(isnan(circ)) = 0;
avg_ecc = mean(ecc,2);
avg_circ = mean(circ,2);
SI_area = mask_area * 0.325*0.325;
SI_perim = mask_perimeter * 0.325;
cell_shape = SI_perim./sqrt(SI_area);
cell_shape(isnan(cell_shape)) = 0;
avg_shape = mean(cell_shape,2);

v(isnan(v)) = 0;
v = v * 0.325;
avg_v = mean(v,2);
idx = find(avg_shape);
c_idx = find(avg_circ);
e_idx = find(avg_ecc);
fitObject = fit(avg_v(idx),avg_shape(idx),'poly1');
figure;plot(fitObject,avg_v(idx),avg_shape(idx), '.k')
xlabel({'Average Single Cell Volocity', '(um/s)'})
ylabel({'Average Cell Shape Index', '(um)'})

figure
plot(avg_v(c_idx),avg_circ(c_idx),'.k')
xlabel({'Average Single Cell Volocity', '(um/s)'})
ylabel({'Average Eccentricity'})

figure
plot(avg_v(e_idx),avg_ecc(e_idx), '.k')
xlabel({'Average Single Cell Volocity', '(um/s)'})
ylabel({'Average Circularity'})
x = avg_v;                                  % Create vxg
y = avg_shape;                    % Create vyg
linearCoefficients = polyfit(x, y, 1);          % Coefficients
yfit = polyval(linearCoefficients, x);          % Estimated  Regression Line
SStot = sum((y-mean(y)).^2);                    % Total Sum-Of-Squares
SSres = sum((y-yfit).^2);                       % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;