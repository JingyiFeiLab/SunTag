%Output structure has one entry for each mRNA. Each entry contains the
%frames, location, and corresponding GFP intensity of that mRNA. After
%running this, you can run "suntag_scatter.m" without changing anything to
%generate images of mRNA spots over GFP images. You can also run
%"sunTag_Intensity.m" without any changes to generate the intensity graphs
%for each mRNA

%Set filepath to folder containing .tif stack
gfp_image_filepath = '/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/IntensityCal71118/5_500';


% Direct to .txt file with mRNA tracking results
gfp_file = '/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/IntensityCal71118/5_500ms.txt';

x_col = 2;
y_col = 1;
z_col = 3;


slices = 1; % How many slices in image
pixelscaling = 130; % Coordinates in nm
window = 1; % Size of mRNA tracking window for GFP projection


particles = textread(gfp_file);

x=particles(:,x_col)/pixelscaling;                          % Select the column containing X coordinates in Pixel form
y=particles(:,y_col)/pixelscaling;                          % Select the column containing Y coordinates in Pixel form
z=particles(:,z_col)/pixelscaling;
spot_frame=particles(:,5);                 % Frame of mRNA Spot
                    % mRNA Spot ID

stack_gfp = imFormat_wholeImage(gfp_image_filepath,slices);
[X,Y] = size(stack_gfp(:,:,1));

x_coord = int32(x);
y_coord = int32(y);
x_coord(x_coord<1) = 1;
y_coord(y_coord<1) = 1;
x_coord(x_coord>=X) = X;
y_coord(y_coord>=Y) = Y;

num_mRNAs = length(particles(:,5));

field1 = 'GFP'; % All Objects, single and multi, labeled
field2 = 'GFP_Center';
field3 = 'Frames';
field4 = 'GFP_Intensity';
field5 = 'GFP_Background';
field6 = 'GFP_noBackSub';
field7 = 'mean_GFP_Pixel';
field8 = 'mean_GFP_Pixel_noBackSub';



GFP_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [],field7, [],field8, []);

for i = 1:num_mRNAs
    GFP_struct(i).GFP = i;
    GFP_struct(i).GFP_Center = [x(i),y(i)];
    frame = particles(i,5)+1;
    GFP_struct(i).Frames = particles(i,5)+1;
    gfp_int_array = zeros(window,window);
    gfp_background_array = zeros(window+2,window+2);
    final_gfp_int_array = [];
    final_gfp_int_array_noBack =[];
    final_gfp_background = [];
  
    final_gfp_pixel = [];
    
    final_gfp_pixel_noBack = [];
    
    temp_x = x_coord(i);
    temp_y = y_coord(i);
    
    center_x_array = temp_x-floor(window/2):temp_x+floor(window/2);
    center_y_array = temp_y-floor(window/2):temp_y+floor(window/2);
    background_x_gfp = temp_x-floor(window/2)-1:temp_x+floor(window/2)+1;
    background_y_gfp = temp_y-floor(window/2)-1:temp_y+floor(window/2)+1;
    background_x_gfp(background_x_gfp<1) = 1;background_y_gfp(background_y_gfp<1) = 1;background_x_gfp(background_x_gfp>X) = X;background_y_gfp(background_y_gfp>Y) = Y;
    center_x_array(center_x_array<1) = 1;center_y_array(center_y_array<1) = 1;center_x_array(center_x_array>X) = X;center_y_array(center_y_array>Y) = Y;
    
    gfp_background_int = zeros(window+2,window+2);
    for ksb = 1:window+2
        for lsb = 1:window
            if sum(ksb==2:window+1) > 0 && sum(lsb==2:window+1) > 0
                continue
            end
            gfp_background_int(ksb,lsb) = stack_gfp(background_y_gfp(ksb),background_x_gfp(lsb),frame);
        end
    end
    
    gfp_background_int = gfp_background_int(:);
    gfp_background_int(gfp_background_int==0) = [];
    gfp_background = median(gfp_background_int);
    
    center_x_array(center_x_array<1) = 1;center_y_array(center_y_array<1) = 1;center_x_array(center_x_array>X) = X;center_y_array(center_y_array>Y) = Y;
    gfp_int = zeros(window,window);
    gfp_int_noBack = zeros(window,window);
    
    for k = 1:window
        for l = 1:window
            q = q + 1;
            gfp_int(k,l) = stack_gfp(center_y_array(k),center_x_array(l),frame) - gfp_background;
            gfp_int_noBack(k,l) = stack_gfp(center_y_array(k),center_x_array(l),frame);
        end
    end
    final_gfp_int_array = sum(gfp_int(:));
    
    final_gfp_int_array_noBack = sum(gfp_int_noBack(:));
    
    final_gfp_pixel = mean(gfp_int(:));
    
    final_gfp_pixel_noBack = mean(gfp_int_noBack(:));
    
    final_gfp_background = gfp_background;
    
    
    GFP_struct(i).GFP_Intensity = final_gfp_int_array;
    
    GFP_struct(i).GFP_noBackSub = final_gfp_int_array_noBack;
    
    GFP_struct(i).GFP_Background = final_gfp_background;
    
    GFP_struct(i).mean_GFP_Pixel = final_gfp_pixel;
    
    GFP_struct(i).mean_GFP_Pixel_noBackSub = final_gfp_pixel_noBack;
    
    
end
                
    
    
    












