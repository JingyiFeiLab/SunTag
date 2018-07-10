
% Change all the names of these directories if you want to save them
% automatically. Alternatively comment them out, as well as lines 12-20,
% and 59-66 if you don't want to change the folder names and save the
% figures automatically

parentDir = 'C:\Users\Jingyi fe\Documents\Data\RNA modification\suntag imaging\07062018_tracking\mut_003002\mut3002\'; %This and line below define
newFile = 'intensity';             % where you want to save the images. Create a new folder
graph_title = 'GFP Intensity for mRNA ';
sample = 'mRNA\'; % Will save images in defined folder, sub-folder called "/mRNA"

background_sub = 1; % 1 = plot background subtacted pixel values. 0 = Plot raw values

compareDir = strcat([parentDir,newFile]);
sampleDir = strcat([compareDir,'\',sample]);
if exist(compareDir,'dir')~=7
    mkdir(parentDir, newFile)
end

if exist(sampleDir,'dir')~=7
    mkdir(compareDir, sample)
end



si = 1;
num_mRNAs = 110;

for i = 1:num_mRNAs
    
    subplot(3,3,si)
    frames = mRNA_struct(i).Frames;
    if background_sub == 1
        gfp = mRNA_struct(i).mean_GFP_Pixel;
        mrna = mRNA_struct(i).mean_mRNA_Pixel;
    elseif background_sub == 0
        gfp = mRNA_struct(i).mean_GFP_Pixel_noBackSub;
        mrna = mRNA_struct(i).mean_mRNA_Pixel_noBackSub;
    end
    
    ylabels{1}='mRNA Pixel Intensity';
    ylabels{2}='GFP Pixel Intensity';
    [ax,hlines1,hlines2] = plotyy(frames,mrna,frames,gfp);
    cfig = get(gcf,'color');
    pos = [0.1  0.1  0.7  0.8];
    offset = pos(3)/5.5;
    set(hlines1,'LineWidth',4);
    set(hlines2,'LineWidth',4)
    set(hlines1,'Color','r')
    set(hlines2,'Color','b')
    if si == 1 || si == 4 || si == 7
        set(get(ax(1),'ylabel'),'string',ylabels{1},'FontWeight','bold')
    end
    if si == 3 || si == 6|| si == 9
        set(get(ax(2),'ylabel'),'string',ylabels{2},'FontWeight','bold')
    end
    set(ax,{'ycolor'},{'r';'b'})
    set(ax,'XLim',[frames(1),frames(length(frames))]);
    %set(ax(1),'YLim',[0 1E5])
    %set(ax(2),'YLim',[0 18E3])
    
    
    
    title(strcat([graph_title, num2str(i)]),'FontSize',14)
    xlabel('Frame','FontSize',14)
    set(ax(1),'FontSize',12);
    set(ax(2),'FontSize',12);
    si = si+1;
    
    if si == 10
        si = 1;
        file1 = strcat([sampleDir,'mRNA',num2str(i)]);
        set(gcf,'position',[315   359   868   623])
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([sampleDir,'mRNA',num2str(i),'.fig']);
        savefig(gcf,file1_fig)
        close all
    end
    
end

%save final figure
file1 = strcat([sampleDir,'mRNA',num2str(i)]);
set(gcf,'position',[315   359   868   623])
set(gcf,'PaperPositionMode','auto')
print(file1,'-painters','-depsc','-r0')
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
file1_fig = strcat([sampleDir,'mRNA',num2str(i),'.fig']);
savefig(gcf,file1_fig)
