parentDir = '/Users/reyer/Data/SunTag_Tracking/'; %This and line below define
newFile = 'IntensityCal71118';             % where you want to save the images. Create a new folder
graph_title = 'mRNA spots over GFP Image, Frame';
sample = 'sample1_5by5';

window = 5;

frames_of_interest = [1];

compareDir = strcat([parentDir,newFile]);
sampleDir = strcat([compareDir,'/',sample]);
if exist(compareDir,'dir')~=7
    mkdir(parentDir, newFile)
end

if exist(sampleDir,'dir')~=7
    mkdir(compareDir, sample)
end

for f = frames_of_interest
    frame = f;
    centers = [];
    
    for i = 1:num_mRNAs
        frames = GFP_struct(i).Frames;
        
        if frames == frame
            centers = [centers; GFP_struct(i).GFP_Center(1),GFP_struct(i).GFP_Center(2)];
        end
        
    end
    
    figure(1);imshow(mat2gray(stack_gfp(:,:,frame)));hold on; scatter(centers(:,1),centers(:,2),50)
    
    
    title(strcat([graph_title,' ',num2str(frame)]),'FontSize',24)
    file1 = strcat([sampleDir,'/frame',num2str(frame)]);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([sampleDir,'/frame',num2str(frame),'.fig']);
    savefig(gcf,file1_fig)
    
    close all
end
