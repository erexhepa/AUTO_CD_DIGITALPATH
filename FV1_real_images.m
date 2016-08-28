clear all
close all

frame = 1;

%% Read the TIF file

name{3} = 'CellBorder1.tif';
name{2} = 'Centriole1.tif';
name{1} = 'Neuron1.tif';
name{16}='Cell_border1.tif';

%  addpath(genpath([cd '/prtools']));
for kid=[1];
    mkdir([cd filesep 'Results' filesep strrep(name{kid},'.tif','')])
    nametex=[cd filesep 'Results' filesep strrep(name{kid},'.tif','') filesep strrep(name{kid},'.tif','_')];
    
    fname=name{kid};
    
    info = imfinfo([cd filesep fname]);
     num_images = numel(info);
    
    Img1=[];
    for k = 1:num_images
        I = imread([cd filesep fname], k);
        Img1(:,:,k)=I;
    end
    %Img1=Img1(1:150,1:150,:);
%
    
    [Img11,zval11]=max(Img1,[],3);
    
    Img11=uint16(65535*mat2gray(Img11));
    
    imwrite(uint16(65535*mat2gray(Img11)),[nametex 'channel1max.png']);
    
    
    IM4=Img1;
    Img2=Img1;
    
    M = [-1 2 -1];
    fname=strrep(fname,'.tif','_');
    
    [zprojf1,qzr2,classmap,idmaxini,cost,WW,sftz2,decay]=Main_SME_method(Img1); 
    
    zprojf1=uint16(65535*(mat2gray(zprojf1)));
    composite_image=zprojf1;
    
    %                 f=strcat(nametex,'composite',num2str(2),'.png');
    %                 imwrite(composite_image,f);
    
    %                      f=strcat(nametex,'map',num2str(2),'.png');
    %                 imwrite(uint16(65535*mat2gray(classmap)),f);
    
    f=strcat(nametex,'zmap2',num2str(2),'.png');
    imwrite(uint8(qzr2*(256/num_images)),f);
    
    zmap=round(qzr2);
    zmap(zmap>k)=k;
    zmap(zmap<1)=1;
    zprojf1=FV1_make_projection_from_layer(Img1,zmap,0,0);
    % zprojf1=FV1_make_projection_from_layer(Img1,zmap,5,5);
    % zprojf3=Img33;
    
    % h = fspecial('gaussian',[5 5],2);
    %                                        zprojf2 = imfilter(zprojf2,h,'replicate');
    %                                        zprojf1 = imfilter(zprojf1,h,'replicate');
    
    
    imwrite(uint16(65535*mat2gray(zprojf1)),[nametex 'channel1proj.png']);
    
    
    
    % f=[nametex 'iteration.tif'];
    %
    % imwrite(movmat(:,:,1)/num_images,f);
    %
    % for kin=2:size(movmat,3)
    %
    %     imwrite(movmat(:,:,kin)/num_images,f,'writemode', 'append');
    %     kin
    % end
    
    
    
    % imwrite(uint16(65535*mat2gray(zprojf2)),[nametex 'channel2proj.tif']);
    
    % IMM=[];
    % IMM(:,:,2)=uint16(65535*mat2gray(zprojf1));
    % IMM(:,:,1)=uint16(65535*mat2gray(zprojf2));
    % IMM(:,:,3)=Img33;
    %
    % imwrite(uint16(65535*mat2gray(IMM)),[nametex 'Maxporjcol2.png']);
    
    
    %# preallocate
    
    % figure, set(gcf, 'Color','white')
    % % Z = peaks; surf(Z);
    % axis tight
    % set(gca, 'nextplot','replacechildren', 'Visible','off');
    %
    % nFrames = size(movmat,3);
    % mov(1:nFrames) = struct('cdata',[], 'colormap',[]);
    
    
    
    % for k = 1:size(movmat,3)
    %     G=movmat(:,:,k);
    %
    % colormap(jet)
    % % subplot(1,2,1)
    % imagesc(G);
    % axis tight
    %   caxis manual
    %                         caxis([1 num_images]);
    %                         colorbar
    % % C = colormap;
    % % L = size(C,1);
    % % Gs = round(interp1(linspace(1,num_images,L),1:L,G));
    % % H = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
    % % subplot(1,2,2)
    % % image(H)  % Does this image match the other one?
    %
    % %                     mov(1,k).cdata =H;
    %
    %
    %
    %                     %# figure
    %
    %
    %
    %
    %    mov(k) = getframe(gca);
    %
    % end
    % close(gcf)
    % %# save as AVI file, and open it using system video player
    % movie2avi(mov, 'myPeaks2.avi','fps',24);
    % % winopen('myPeaks1.avi')
    %
    close all
end

