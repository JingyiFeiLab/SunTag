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
slices = 31;
pixelscaling = 1;
window = 5; % Size of mRNA tracking window for GFP projection


particles = textread(mRNA_file);

x=particles(:,x_col)/pixelscaling;                          % Select the column containing X coordinates in Pixel form
y=particles(:,y_col)/pixelscaling;                          % Select the column containing Y coordinates in Pixel form
spot_frame=particles(:,4);                 % Frame of mRNA Spot
SpotID=particles(:,5);                     % mRNA Spot ID

stack_gfp = imFormat_wholeImage(gfp_image_filepath,slices);
stack_mRNA = imFormat_wholeImage(mRNA_image_filepath,slices);
[X,Y] = size(stack_gfp(:,:,1));


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


mRNA_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [],field7, [],field8, []);

for i = 1:num_mRNAs
    mRNA_struct(i).mRNA = i;
    temp_mRNA = particles(particles(:,5)==i,:);
    mRNA_struct(i).mRNA_Center = [temp_mRNA(:,x_col),temp_mRNA(:,y_col)];
    frames = temp_mRNA(:,4)+1;
    mRNA_struct(i).Frames = frames;
    gs = 10-window;
    gfp_int_array = zeros(gs^2,length(frames));
    gfp_center = zeros(length(frames),2);
    final_gfp_int_array = zeros(1,length(frames));
    final_mRNA_int_array = zeros(1,length(frames));
    final_gfp_background = zeros(1,length(frames));
    final_mRNA_background = zeros(1,length(frames));
    for j = 1:length(frames)
        temp_x = int32(temp_mRNA(j,x_col));
        temp_y = int32(temp_mRNA(j,y_col));
        center_x_mRNA = temp_x-1:temp_x+1;
        center_y_mRNA = temp_y-1:temp_y+1;
        background_x_mRNA = temp_x-3:temp_x+3;
        background_y_mRNA = temp_y-3:temp_y+3;
        background_x_mRNA(background_x_mRNA<1) = 1;background_y_mRNA(background_y_mRNA<1) = 1;background_x_mRNA(background_x_mRNA>X) = X;background_y_mRNA(background_y_mRNA>Y) = Y;
        
        mRNA_background_int = zeros(7,7);
        for msb = 1:7
            for nsb = 1:7
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
        
        qs = 0;
        for ms = 1:3
            for ns = 1:3
                qs = qs + 1;
                mRNA_int(ms,ns) = stack_mRNA(center_y_mRNA(ms),center_x_mRNA(ns),frames(j)) - mRNA_background;
            end
        end
        
                
        
        center_x_array = temp_x-floor(gs/2):temp_x+floor(gs/2);
        center_y_array = temp_y-floor(gs/2):temp_y+floor(gs/2);
        background_x_gfp = temp_x-5:temp_x+5;
        background_y_gfp = temp_y-5:temp_y+5;
        background_x_gfp(background_x_gfp<1) = 1;background_y_gfp(background_y_gfp<1) = 1;background_x_gfp(background_x_gfp>X) = X;background_y_gfp(background_y_gfp>Y) = Y;
        
        gfp_background_int = zeros(11,11);
        for ksb = 1:11
            for lsb = 1:11
                if sum(ksb==2:10) > 0 && sum(lsb==2:10) > 0 
                    continue
                end
                gfp_background_int(ksb,lsb) = stack_gfp(background_y_gfp(ksb),background_x_gfp(lsb),frames(j));
            end
        end
        
        gfp_background_int = gfp_background_int(:);
        gfp_background_int(gfp_background_int==0) = [];
        gfp_background = median(gfp_background_int);
        
        
        q = 0;
        for m = 1:gs
            for n = 1:gs
                q = q + 1;
                center_x = center_x_array(m);
                center_y = center_y_array(n);
                temp_x_array = center_x-floor(window/2):center_x+floor(window/2);
                temp_y_array = center_y-floor(window/2):center_y+floor(window/2);
                temp_x_array(temp_x_array<1)=1;temp_y_array(temp_y_array<1)=1;temp_x_array(temp_x_array>X)=X;temp_y_array(temp_y_array>Y)=Y;
                gfp_int = zeros(5,5);
                for k = 1:window
                    for l = 1:window
                        gfp_int(k,l) = stack_gfp(temp_y_array(k),temp_x_array(l),frames(j)) - gfp_background;
                    end
                end
                gfp_int_array(q,j) = sum(gfp_int(:));
            end
        end
        [final_gfp_int_array(j),gfp_ind] = max(gfp_int_array(:,j));
        [gfp_ind_x,gfp_ind_y] = ind2sub(size(zeros(5,5)),gfp_ind);
        final_mRNA_int_array(j) = sum(mRNA_int(:));
        final_gfp_background(j) = gfp_background;
        final_mRNA_background(j) = mRNA_background;
        gfp_center(j,:) = [center_x_array(gfp_ind_x),center_y_array(gfp_ind_y)];
    end
    final_gfp_int_array(final_gfp_int_array<0) = 0;
    final_mRNA_int_array(final_mRNA_int_array<0) = 0;
    mRNA_struct(i).gfp_Center = gfp_center;
    mRNA_struct(i).GFP_Intensity = final_gfp_int_array';
    mRNA_struct(i).mRNA_Intensity = final_mRNA_int_array';
    mRNA_struct(i).GFP_Background = final_gfp_background';
    mRNA_struct(i).mRNA_Background = final_mRNA_background';
    
end
                
    
    
    












