%Output structure has one entry for each mRNA. Each entry contains the
%frames, location, and corresponding GFP intensity of that mRNA. After
%running this, you can run "suntag_scatter.m" without changing anything to
%generate images of mRNA spots over GFP images. You can also run
%"sunTag_Intensity.m" without any changes to generate the intensity graphs
%for each mRNA

%Set filepath to folder containing .tif stack
gfp_image_filepath = '/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/suntag8418/c2';
mRNA_image_filepath = '/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/suntag8418/c1';


% Direct to .txt file with mRNA tracking results
mRNA_file = '/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/suntag8418/sample8418.txt';

x_col = 2;
y_col = 1;
z_col = 3;

% Don't Change any of these
channels = 1;
ref_channel = 1;
ref_slice = 16;
dim = 2;
slices2D = 15;
pixelscaling = 1;
window = 1; % Size of mRNA tracking window for GFP projection




particles = textread(mRNA_file);

x=particles(:,x_col)/pixelscaling;                          % Select the column containing X coordinates in Pixel form
y=particles(:,y_col)/pixelscaling;                          % Select the column containing Y coordinates in Pixel form
spot_frame=particles(:,4);                 % Frame of mRNA Spot
SpotID=particles(:,5);                     % mRNA Spot ID

[~,~,stack_gfp,~,~,~,~,slices] = imFormat(gfp_image_filepath,ref_channel,dim,ref_slice,slices2D,channels);
[~,~,stack_mRNA,~,~,~,~,slices] = imFormat(mRNA_image_filepath,ref_channel,dim,ref_slice,slices2D,channels);
[X,Y] = size(stack_gfp(:,:,ref_slice));

x_coord = int32(x);
y_coord = int32(y);
x_coord(x_coord<1) = 1;
y_coord(y_coord<1) = 1;
x_coord(x_coord>=X) = X;
y_coord(y_coord>=Y) = Y;

num_mRNAs = max(SpotID);

field1 = 'mRNA'; % All Objects, single and multi, labeled
field2 = 'mRNA_Center';
field3 = 'gfp_Center';
field4 = 'Frames';
field5 = 'GFP_Intensity';
field6 = 'mRNA_Intensity';
field7 = 'GFP_Background';
field8 = 'mRNA_Background';
field9 = 'GFP_noBackSub';
field10 = 'mRNA_noBackSub';
field11 = 'mean_GFP_Pixel';
field12 = 'mean_mRNA_Pixel';
field13 = 'mean_GFP_Pixel_noBackSub';
field14 = 'mean_mRNA_Pixel_noBackSub';


mRNA_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [],field7, [],field8, [],field9, [],field10, [],field11, [],field12, [],field13, [],field14, []);

for i = 1:num_mRNAs
    mRNA_struct(i).mRNA = i;
    temp_mRNA = particles(particles(:,5)==i,:);
    mRNA_struct(i).mRNA_Center = [temp_mRNA(:,x_col),temp_mRNA(:,y_col)];
    mRNA_struct(i).gfp_Center = [temp_mRNA(:,x_col),temp_mRNA(:,y_col)];
    frames = temp_mRNA(:,4)+1;
    mRNA_struct(i).Frames = frames;
    gfp_int_array = zeros(window,window);
    gfp_background_array = zeros(window+2,window+2);
    final_gfp_int_array = zeros(1,length(frames));
    final_mRNA_int_array = zeros(1,length(frames));
    final_gfp_int_array_noBack = zeros(1,length(frames));
    final_mRNA_int_array_noBack = zeros(1,length(frames));
    final_gfp_background = zeros(1,length(frames));
    final_mRNA_background = zeros(1,length(frames));
    
    final_gfp_pixel = zeros(1,length(frames));
    final_mRNA_pixel = zeros(1,length(frames));
    final_gfp_pixel_noBack = zeros(1,length(frames));
    final_mRNA_pixel_noBack = zeros(1,length(frames));
    
    for j = 1:length(frames)
        temp_x = int32(temp_mRNA(j,x_col));
        temp_y = int32(temp_mRNA(j,y_col));
        center_x_mRNA = temp_x-1:temp_x+1;
        center_y_mRNA = temp_y-1:temp_y+1;
        background_x_mRNA = temp_x-3:temp_x+3;
        background_y_mRNA = temp_y-3:temp_y+3;
        background_x_mRNA(background_x_mRNA<1) = 1;background_y_mRNA(background_y_mRNA<1) = 1;background_x_mRNA(background_x_mRNA>X) = X;background_y_mRNA(background_y_mRNA>Y) = Y;
        
        mRNA_background_int = zeros(5,5);
        for msb = 1:5
            for nsb = 1:5
                if sum(background_x_mRNA(msb)==center_x_mRNA)>0 && sum(background_y_mRNA(nsb)==center_y_mRNA)>0
                    continue
                end
                mRNA_background_int(msb,nsb) = stack_mRNA(background_y_mRNA(msb),background_x_mRNA(nsb),frames(j));
            end
        end
        
        mRNA_background_int = mRNA_background_int(:);
        mRNA_background_int(mRNA_background_int==0) = [];
        mRNA_background = median(mRNA_background_int);
        
        center_x_mRNA(center_x_mRNA<1) = 1;center_y_mRNA(center_y_mRNA<1) = 1;center_x_mRNA(center_x_mRNA>X) = X;center_y_mRNA(center_y_mRNA>Y) = Y;
        mRNA_int = zeros(3,3);
        mRNA_int_noBack = zeros(3,3);
        
        for ms = 1:3
            for ns = 1:3
                
                mRNA_int(ms,ns) = stack_mRNA(center_y_mRNA(ms),center_x_mRNA(ns),frames(j)) - mRNA_background;
                mRNA_int_noBack(ms,ns) = stack_mRNA(center_y_mRNA(ms),center_x_mRNA(ns),frames(j));
            end
        end
        
                
        
        center_x_array = temp_x-floor(window/2):temp_x+floor(window/2);
        center_y_array = temp_y-floor(window/2):temp_y+floor(window/2);
        background_x_gfp = temp_x-floor(window/2)-1:temp_x+floor(window/2)+1;
        background_y_gfp = temp_y-floor(window/2)-1:temp_y+floor(window/2)+1;
        background_x_gfp(background_x_gfp<1) = 1;background_y_gfp(background_y_gfp<1) = 1;background_x_gfp(background_x_gfp>X) = X;background_y_gfp(background_y_gfp>Y) = Y;
        
        gfp_background_int = zeros(window+2,window+2);
        for ksb = 1:window+2
            for lsb = 1:window
                if sum(ksb==2:window+1) > 0 && sum(lsb==2:window+1) > 0
                    continue
                end
                gfp_background_int(ksb,lsb) = stack_gfp(background_y_gfp(ksb),background_x_gfp(lsb),frames(j));
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
                gfp_int(k,l) = stack_gfp(center_y_array(k),center_x_array(l),frames(j)) - gfp_background;
                gfp_int_noBack(k,l) = stack_gfp(center_y_array(k),center_x_array(l),frames(j));
            end
        end
        
        final_gfp_int_array(j) = sum(gfp_int(:));
        final_mRNA_int_array(j) = sum(mRNA_int(:));
        final_gfp_int_array_noBack(j) = sum(gfp_int_noBack(:));
        final_mRNA_int_array_noBack(j) = sum(mRNA_int_noBack(:));
        final_gfp_pixel(j) = mean(gfp_int(:));
        final_mRNA_pixel(j) = mean(mRNA_int(:));
        final_gfp_pixel_noBack(j) = mean(gfp_int_noBack(:));
        final_mRNA_pixel_noBack(j) = mean(mRNA_int_noBack(:));
        
        final_gfp_background(j) = gfp_background;
        final_mRNA_background(j) = mRNA_background;
    end
%     final_gfp_int_array(final_gfp_int_array<0) = 0;
%     final_mRNA_int_array(final_mRNA_int_array<0) = 0;
    mRNA_struct(i).GFP_Intensity = final_gfp_int_array';
    mRNA_struct(i).mRNA_Intensity = final_mRNA_int_array';
    mRNA_struct(i).GFP_noBackSub = final_gfp_int_array_noBack';
    mRNA_struct(i).mRNA_noBackSub = final_mRNA_int_array_noBack';
    mRNA_struct(i).GFP_Background = final_gfp_background';
    mRNA_struct(i).mRNA_Background = final_mRNA_background';
    mRNA_struct(i).mean_GFP_Pixel = final_gfp_pixel';
    mRNA_struct(i).mean_mRNA_Pixel = final_mRNA_pixel';
    mRNA_struct(i).mean_GFP_Pixel_noBackSub = final_gfp_pixel_noBack';
    mRNA_struct(i).mean_mRNA_Pixel_noBackSub = final_mRNA_pixel_noBack';
    
end
                
    
    
    












