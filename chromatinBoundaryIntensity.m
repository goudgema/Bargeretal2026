clear

%Written originally April 2024, updated April 2026 for submission.
%This takes in a a time course 
%image with signal on the nuclear envelope (int_name)
%and a premasked (generated in ImageJ, curv_name) chromatin mass time course
%and outputs the mean intensity per image and mean variance per image in
%the column variables mean_intensities and variance_intensities

%file declerations
curv_name = "TOFA24m.tif";
int_name= "TOFA24.tif";
number_images = max(size(imfinfo(curv_name)));

%empty decleration
mean_intensities = zeros(number_images,1);
variance_intensities = zeros(number_images,1);

for im_num = 1:number_images

    %image declerations for each iteration
    curv_image = double(imread(curv_name, im_num));
    [B, L] = bwboundaries(curv_image);
    curv_image = zeros(size(curv_image));

    for k = 1:length(B)
    boundary = B{k}; % Nx2 array [row, col]
    
    % Convert subscripts to linear indices
    idx = sub2ind(size(curv_image), boundary(:,1), boundary(:,2));
    
    % Set those pixels to 1
    curv_image(idx) = 1;
    end
    int_image = double(imread(int_name, im_num));


    %inefficient but kept from how the old code was constructed based on
    %boundaries initially.
    binary_curv = curv_image;
    binary_curv(binary_curv ~= 0) = 1;
    filled_curv = imfill(binary_curv, "holes");
    
    %distance transform, change the distnace > 2 for desired pixel width
    distance = bwdist(filled_curv);
    temp_int = int_image;
    temp_int(distance > 2 | distance == 0) = 0;

    %iterates through the pixels that comply within the distance from
    %boundary
    [nonZeroRows, nonZeroColumns] = find(temp_int);
    number_of_pixels = max(size(nonZeroColumns));
    number_in_boundary = max(size(nonzeros(binary_curv)));
    intensity = zeros(number_of_pixels+number_in_boundary, 1);
    
    for k = 1:number_of_pixels

        %makke a new distance transform for that pixel
        blank_distance = zeros(size(binary_curv));
        blank_distance(nonZeroRows(k), nonZeroColumns(k)) = 1;
        distance_transform = bwdist(blank_distance);
        distance_transform(binary_curv == 0) = NaN;
        min_dist = min(min(distance_transform, [], "omitnan"), [], "omitnan");
        [closest_y, closest_x] = find(distance_transform == min_dist);
        intensity(k) = int_image(nonZeroRows(k), nonZeroColumns(k));
    end
    

    %redeclare initialize to now loop through all the points on the curve
    [nonZeroRows, nonZeroColumns] = find(binary_curv);
    
    for m = (number_of_pixels + 1):(number_of_pixels + number_in_boundary)
        corrected_index = m-number_of_pixels;
        intensity(m) = int_image(nonZeroRows(corrected_index), nonZeroColumns(corrected_index));
    end
    
    %background adjdust, manual done for shared-scope 
    intensity = intensity-120;
    mean_intensities(im_num) = mean(intensity);
    variance_intensities(im_num) = var(intensity);
end





