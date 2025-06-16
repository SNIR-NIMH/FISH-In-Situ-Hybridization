function [img2_reg,opts]=FISH_registration(img1,img2,outputdir,varargin)
% 
% function [moving_reg,opts] = FISH_registration(fixed,moving,outputdir,varargin)
% 
% FIXED       Source/target/fixed image, or a folder containing multiple 2D tifs
% MOVING      Moving image, to be registered to the target fixed image. It can
%             also be a folder containing same number of 2D tifs as FIXED, where
%             each of them will be registered to the corresponding tif in FIXED
% OUTPUTDIR   Output directory where all output files will be written
% 
% Optional name/value parameters
% CHUNKSIZE   (Optional) Chunksize, e.g. 10x10, a string separated by x. If not
%             mentioned, it will automatically be computed.
% DSFACTOR    (Optional) Initial downsampling factor to try the coarse
%             registration. Too high downsampling factor cause unstable
%             registration. Enter a range in Matlab notation in increasing
%             order, such as 50:-2:10. Default is 31:-2:3.
% Overlap     Overlap factor between chunks. Default is 0.05 (5% overlap). Enter
%             a number between 0 and 1.
% MaxTranslationPixels
%             (Optional) Maximum translation in pixels that will be tolerated
%             during the chunk-wise registration. If not mentioned, default is 
%             2x overlap, i.e., 10% of the chunk size. Increase this tolerance
%             if there is significant translation in the image.
% WriteOutput (Optional) Either 1 or 0. Default 1, indicating the registered
%             image and transformation matrix will be written to disk as a tif
%             file and a mat file, respectively. 0 indicates only their value
%             will be returned (useful in interactive mode).
% InitReg     (Optional) Initial registration type to estimate the proper
%             downsampling factor. Default is rigid. Options are rigid or
%             translation.
% FinalReg    (Optional) Registration type to use to register chunks. Default is
%             similarity (i.e. affine). Options are similarity (affine), rigid, 
%             or translation.
% Modality    (Optional) Either monomodal or multimodal. Default monomodal.
% NumCPU      (Optional) When both FIXED and MOVING are folders, they can be
%             processed in parallel using NumCPU processes. Default 12.
% Pyramidlevel  (Optional) Number of pyramid levels to use for registraton.
%             Default is 2. Use higher (4) if the images are large. Using small
%             (1) may incur in bad registration. This is used only when modality
%             is multimodal.
% Mask        (Optional) A mask image to introduce better registration  when
%             there are a lot of background noise
% Quantile    (Optional) Default 0.99. An upper quantile to clamp the outlier
%             intensities. It is useful if there are a lot of outliers that can
%             cause registration problems.


% NOTE: For multimodal registration, imregtform works better than imregcorr. I will
% use that for modality = multimodal

warning off;
% tic
ds=20:-2:3; % default.
chunks=[-1 -1]; % default is optimal chunking based on registration
tError=-1; % default 2*ovl (defined later)
ovl=0.05; % default 5% overlap
write=1;  % default write the output image and transformation matrix
initregtype='rigid';
finalregtype='similarity';
modal='monomodal';
numcpu=12;
img2_reg=[];
opts=[];
pyramid=2;
mask=[];
quant=0.99;
for i = 1:2:length(varargin)
    s=lower(varargin{i});
    if ~strcmpi(s,'chunksize') && ~strcmpi(s,'dsfactor') ...
            && ~strcmpi(s,'maxtranslationpixels') && ~strcmpi(s,'overlap') && ...
            ~strcmpi(s,'writeoutput') && ~strcmpi(s,'initreg')  && ~strcmpi(s,'modality') ...
             && ~strcmpi(s,'finalreg')  && ~strcmpi(s,'numcpu') ...
             && ~strcmpi(s,'mask')   && ~strcmpi(s,'pyramidlevel') && ~strcmpi(s,'quantile' ) 
        fprintf('ERROR: Argument \"%s\" is not recognized.\n',varargin{i});
        return;
    end
    if strcmpi(s,'chunksize') % if chunksize is not mentioned it will automatically be created
        if ~ischar(varargin{i+1})
            fprintf('ERROR: Chunksize must be a string separated by x, e.g. 10x10.\n');
            return;
        end
            
        temp=strsplit(varargin{i+1},'x');
        if length(temp)~=2
            fprintf('ERROR: Chunksize must be a string separated by x, e.g. 10x10.\n');
            return;
        end
        chunks=zeros(1,2);
        for t=1:2
            chunks(t)=str2num(temp{t});
        end
    elseif strcmpi(s,'modality')  % downsampling factor to check for robust registration
        modal=varargin{i+1};
        if ~strcmpi(modal,'monomodal') && ~strcmpi(modal,'multimodal')
            fprintf('ERROR: modality must be either monomodal or multimodal.\n');
            return;
        end
        
    elseif strcmpi(s,'dsfactor')  % downsampling factor to check for robust registration
        if isdeployed
            ds=str2num(varargin{i+1});
        else
            ds=varargin{i+1};
        end
    elseif strcmpi(s,'maxtranslationpixels')  % 
        if isdeployed
            tError=str2num(varargin{i+1});
        else
            tError=varargin{i+1};
        end    
    elseif strcmpi(s,'overlap')  % overlap factor between image.
        if isdeployed
            ovl=str2num(varargin{i+1});
        else
            ovl=varargin{i+1};
        end 
        if ovl<=0 || ovl>=1
            fprintf('ERROR: Overlap must be a number between 0 and 1.\n');
            return;
        end
    elseif strcmpi(s,'writeoutput')  % apply the transformation matrix on a new image
        if isdeployed
            write=str2num(varargin{i+1});
        else
            write=varargin{i+1};
        end 
        if write~=0 && write~=1
            fprintf('ERROR: WriteOutput must be either 0 or 1.\n');
            return;
        end
    elseif strcmpi(s,'initreg')  % apply the transformation matrix on a new image
        initregtype=lower(varargin{i+1});
        if ~strcmpi(initregtype,'rigid') && ~strcmpi(initregtype,'translation') ...
                && ~strcmpi(initregtype,'similarity')
            fprintf('ERROR: InitReg must be rigid/translation/similarity.\n');
            return;
        end
    elseif strcmpi(s,'finalreg')  % apply the transformation matrix on a new image
        finalregtype=lower(varargin{i+1});
        if ~strcmpi(finalregtype,'rigid') && ~strcmpi(finalregtype,'translation') ...
                && ~strcmpi(finalregtype,'similarity')
            fprintf('ERROR: FinalReg must be rigid/translation/similarity(affine).\n');
            return;
        end
    elseif strcmpi(s,'numcpu')  % downsampling factor to check for robust registration
        if isdeployed
            numcpu=str2num(varargin{i+1});
        else
            numcpu=varargin{i+1};
        end
    elseif strcmpi(s,'mask')  
        mask=varargin{i+1};
        if ~isfile(mask) 
            fprintf('ERROR: The mask must be a 2D tif file.');
            return;
        end
    elseif strcmpi(s,'pyramidlevel')  % downsampling factor to check for robust registration
        if isdeployed
            pyramid=str2num(varargin{i+1});
        else
            pyramid=varargin{i+1};
        end 
    elseif strcmpi(s,'quantile')  % quantile to clamp image intensities
        if isdeployed
            quant=str2num(varargin{i+1});
        else
            quant=varargin{i+1};
        end
    else
        fprintf('ERROR: Argument \"%s\" is not acceptable.\n',s);
        return;
    end  
    
end



if ~isfolder(outputdir)
    mkdir(outputdir);
end
A=rdir(fullfile(outputdir,'*.tif'));
if ~isempty(A)
    fprintf('WARNING: Output folder contains some tif images. They may be overwritten.\n');
end

args.ds=ds; % default.
args.chunks=chunks; % default is optimal chunking based on registration
args.tError=tError; % default 2*ovl (defined later)
args.ovl=ovl; % default 5% overlap
args.write=write;  % default write the output image and transformation matrix
args.initregtype=initregtype;
args.finalregtype=finalregtype;
args.modal=modal;
args.numcpu=numcpu;
args.outputdir=outputdir;
args.mask=mask;
args.pyramid=pyramid;
args.quant=quant;
disp(args);

if isfile(img1) && isfile(img2)
    [img2_reg,opts]=register_pair_of_images(img1,img2,args);
elseif isfolder(img1) && isfolder(img2)
    filelist1=rdir(fullfile(img1,'*.tif'));
    filelist2=rdir(fullfile(img2,'*.tif'));
    if length(filelist1) ~= length(filelist2)
        fprintf('ERROR: Number of images in FIXED folder (%d) must match with that of MOVING folder (%d).\n',...
            length(filelist1),length(filelist2));
        return;
    end
    if isempty(gcp('nocreate'))
        username=getenv('USER');
        tempdirname=tempname(fullfile('/home',username,'.matlab','local_cluster_jobs','R2022a'));
        mkdir(tempdirname);
        cluster=parallel.cluster.Local();
        cluster.NumWorkers=numcpu;
        cluster.JobStorageLocation=tempdirname;
        fprintf('Temp Job directory = %s\n',tempdirname);
        pl=parpool(cluster);
    else
        pl=[];
    end
    parfor k=1:length(filelist1)
%         img1=imread(filelist1(k).name);  % don't do this, output filename is
%         required
%         img2=imread(filelist2(k).name);
        register_pair_of_images(filelist1(k).name,filelist2(k).name,args);

    end
else
    fprintf('ERROR: FIXED image type (either single 2D tif or a folder) must match with that of MOVING image type.\n');
    return;
end


end


function [img2_reg,opts]=register_pair_of_images(img1,img2,args)
tic
ds=args.ds; % default.
chunks=args.chunks; % default is optimal chunking based on registration
tError=args.tError; % default 2*ovl (defined later)
ovl=args.ovl; % default 5% overlap
write=args.write;  % default write the output image and transformation matrix
initregtype=args.initregtype;
finalregtype=args.finalregtype;
modal=args.modal;
outputdir=args.outputdir;
mask=args.mask;
quant=args.quant;
img2_reg=[];opts=[]; % just initialize


if ischar(img1)
    fprintf('Reading %s\n',img1);
    x=imread(img1);
    x=single(x);
elseif isnumeric(img1)
    x=single(img1);
end



q1=quantile(x(x>0),quant);
x(x>q1)=q1;
if ischar(img2)
    fprintf('Reading %s\n',img2);
    y=imread(img2);
    y=single(y);
elseif isnumeric(img2)
    y=single(img2);
end
origy=y; % keep a copy because of intensity clamping
q2=quantile(y(y>0),quant);
y(y>q2)=q2;
fprintf('Clamping image intensities of fixed and moving images at %.2f and %.2f.\n',q1,q2);

dim1=size(x);
dim2=size(y);
if ~isempty(mask)
    if ischar(mask)
        fprintf('Reading %s\n',mask);
        m=imread(mask);
    elseif isnumeric(img2)
        m=single(mask);
    end
    m=uint8(m>0);
    if size(m,1) ~=dim1(1) || size(m,2) ~=dim1(2)
        fprintf('ERROR: Mask must have the same size the fixed image.\n');
        return;
    end

end

dim=[max(dim1(1),dim2(1)) max(dim1(2),dim2(2))];
x1=zeros(dim,'single');
y1=zeros(dim,'single');
y2=zeros(dim,'single');
% Padding images to have same size
x1(1:dim1(1),1:dim1(2))=x;
y1(1:dim2(1),1:dim2(2))=y;
y2(1:dim2(1),1:dim2(2))=origy;
x=x1;y=y1;clear x1 y1;
origy=y2;clear y2;
if ~isempty(mask)
    m1=zeros(dim,'uint16');
    m1(1:dim1(1),1:dim1(2))=m;
    m=m1;
end
fprintf('Padded image size = %d x %d \n',size(x,1),size(x,2));


% downsample and rigid register with various downsampling factor to get a robust
% rigid registration parameter
T=zeros(length(ds),2);
rotangle=zeros(1,length(ds));
count=1;

fprintf('Downsampling by %d to %d and registering with a %s registration:\n',min(ds),max(ds),initregtype);
[optimizer, metric]=imregconfig(modal);
optimizer.MaximumIterations = 1000;
for dsfactor=ds  % can't use mask here because x and y are not registered
%     dsfactor=ceil(max(dim/siz));
    fprintf('%d,',dsfactor);
    dim2=round(dim(1:2)/dsfactor);
    x1=imresize(x,dim2,'bicubic');
    y1=imresize(y,dim2,'bicubic');
%     if isempty(mask)
%         x1=imresize(x,dim2,'bicubic');
%         y1=imresize(y,dim2,'bicubic');
%     else
%         x1=imresize(x.*m,dim2,'bicubic');
%         y1=imresize(y.*m,dim2,'bicubic');
%     end
%     x1=x(1:dsfactor:end,1:dsfactor:end);
%     y1=y(1:dsfactor:end,1:dsfactor:end);
    if strcmpi(modal,'multimodal')
        tformEstimate=imregtform(y1,x1,initregtype,optimizer,metric,'PyramidLevels',args.pyramid);
    else
        tformEstimate = imregcorr(y1,x1, initregtype,'Window',true);
    end
    tformEstimate.T(3,1:2)=dsfactor*tformEstimate.T(3,1:2);
    T(count,:)=tformEstimate.T(3,1:2);
    E=Rotationmat2Euler(tformEstimate.T(1:2,1:2));
    rotangle(count)=rad2deg(E(3));
%     Rfixed = imref2d(size(x1));
%     y_tf=imwarp(y1,tformEstimate,'OutputView',Rfixed);
    count=count+1;
end
fprintf('\n');
q1=quantile(rotangle,0.25);
q2=quantile(rotangle,0.75);
indx=find(rotangle>=q1 & rotangle<=q2);
rotangle=rotangle(indx);
T=T(indx,:);
ds=ds(indx);

if length(unique(rotangle))==1 % rotation angles are all same, nothing important, use translations
    T1=sum(sqrt(T.^2),2);
    indx1=find(T1==median(T1));
    % For some cases, the median could be outside of the range, choose the one that
    % is closest to the median
    if isempty(indx1)
        d=abs(median(T1)-T1);
        indx1=find(d==min(d));
        
    end
else
    
    indx1=find(rotangle==median(rotangle));
    % For some cases, the median could be outside of the range, choose the one that
    % is closest to the median
    if isempty(indx1)
        d=abs(median(rotangle)-rotangle);
        indx1=find(d==min(d));
        
    end
end

if length(indx1)>1
    fprintf('A robust global transformation matrix found at downsampling factor = %d.\n',ds(indx1(1)));
    indx1=indx1(1);
else
    fprintf('Global transformation matrix found at downsampling factor = %d.\n',ds(indx1));
end
dsfactor=ds(indx1);

if chunks(1)==-1 && chunks(2)==-1
    chunks=[dsfactor dsfactor];
    fprintf('Automatically splitting image into %d x %d tiles for tile-wise registration.\n',chunks(1),chunks(2));
else
   
    fprintf('Splitting image into %d x %d tiles for tile-wise registration.\n',chunks(1),chunks(2));
end

dim2=round(dim(1:2)/dsfactor);
x1=imresize(x,dim2,'bicubic');
y1=imresize(y,dim2,'bicubic');
% if isempty(mask)
%     x1=imresize(x,dim2,'bicubic');
%     y1=imresize(y,dim2,'bicubic');
% else
%     x1=imresize(x.*m,dim2,'bicubic');
%     y1=imresize(y.*m,dim2,'bicubic');
% end
% x1=x(1:dsfactor:end,1:dsfactor:end);
% y1=y(1:dsfactor:end,1:dsfactor:end);
if strcmpi(modal,'multimodal')
    tformInit=imregtform(y1,x1,initregtype,optimizer,metric,'PyramidLevels',args.pyramid);
else
    tformInit = imregcorr(y1,x1,initregtype,'Window',true);
end
tformInit.T(3,1:2)=dsfactor*tformInit.T(3,1:2);
E=Rotationmat2Euler(tformInit.T(1:2,1:2));
fprintf('Transformation matrix (pixels) and z-rotation = %.4f degrees\n',E(3));
disp(tformInit.T);
Rfixed = imref2d(size(x));
y = imwarp(y,tformInit,'OutputView',Rfixed);
origy = imwarp(origy,tformInit,'OutputView',Rfixed);

% Now mask the transformed images
if ~isempty(mask)
    x=x.*single(m);
    y=y.*single(m);
end


N=prod(chunks);
padsize=max(round(ovl*(dim./chunks)));  % a 5% overlap between tiles
% pdim1=round(dim./chunks);  
pdim=floor(dim./chunks)+2*padsize;  % dimension of the overlapping tiles to be registered
tform=zeros([3 3 N],'single');
regY=cell(chunks(1),chunks(2));
% corrvol=cell(chunks(1),chunks(2));
fprintf('Padding size = %d pixels.\n',padsize);
fprintf('Registration chunk size = %d x %d pixels.\n',pdim(1),pdim(2));
d=floor(dim./chunks);
count=1;
if tError==-1
    tError=2*padsize;
end
[optimizer, metric]=imregconfig(modal);
optimizer.MaximumIterations = 800;
if strcmpi(modal,'multimodal')
    metric.NumberOfSpatialSamples=2000;
    metric.NumberOfHistogramBins=256;
end


fprintf('Error tolerance = %d pixels.\n',tError);
for r=1:chunks(1)
    for c=1:chunks(2)
        %%
        tempx=zeros(pdim,'single');
        tempy=zeros(pdim,'single');
        tempyorig=zeros(pdim,'single');
        I1=(r-1)*d(1)+1-padsize;
        I2=r*d(1)+padsize;
        J1=(c-1)*d(2)+1-padsize;
        J2=c*d(2)+padsize;
        if I1<1
            I11=-I1; % index of tempx
            I1=1;    % index of x
            
        else
            I11=1;
        end
        if J1<0
            J11=-J1; % index of tempy
            J1=1;    % index of y           
        else
            J11=1;
        end
        if I2>dim(1)
            I2=dim(1);        
        end
        if J2>dim(2)
            J2=dim(2);        
        end
        a=x(I1:I2,J1:J2);
        adim=size(a);
        tempx(I11:I11+adim(1)-1,J11:J11+adim(2)-1)=single(a);
        
        a=y(I1:I2,J1:J2);        
        tempy(I11:I11+adim(1)-1,J11:J11+adim(2)-1)=single(a);
        a=origy(I1:I2,J1:J2);        
        tempyorig(I11:I11+adim(1)-1,J11:J11+adim(2)-1)=single(a);

        
        % Use imregcorr instead of imregtform for monomodal images
        if strcmpi(modal,'monomodal')
            if sum(tempx(:))>0 && sum(tempy(:))>0  % ignore blank tiles at the
                tformEstimate = imregcorr(tempy,tempx,finalregtype,'Window',true);
                % translation only is sometimes more robust than rigid, but
                % rigid/similarity works for small tiles
                if abs(tformEstimate.T(3,1))>=tError ||abs(tformEstimate.T(3,2))>=tError
                    fprintf('WARNING: Potential failure, displacement = %d,%d pixels. Retrying with %s (window=false).\n',...
                        round(tformEstimate.T(3,1)),round(tformEstimate.T(3,2)),finalregtype);
                    tformEstimate = imregcorr(tempy,tempx,finalregtype,'Window',false);
                    
                    if abs(tformEstimate.T(3,1))>=tError ||abs(tformEstimate.T(3,2))>=tError
                        fprintf('WARNING: Second potential failure, displacement = %d,%d pixels. Retrying with translation (window=false).\n',...
                            round(tformEstimate.T(3,1)),round(tformEstimate.T(3,2)));
                        tformEstimate = imregcorr(tempy,tempx,'translation','Window',false);
                        
                        if abs(tformEstimate.T(3,1))>=tError ||abs(tformEstimate.T(3,2))>=tError
                            fprintf('WARNING: Third potential failure, displacement = %d,%d pixels. Retrying with imregtform.\n',...
                                round(tformEstimate.T(3,1)),round(tformEstimate.T(3,2)));
                            tformEstimate=imregtform(tempy,tempx,'translation',optimizer,metric,'PyramidLevels',args.pyramid);
                            
                        end
                    end
                    
                end
                flag=1;
            else
                tformEstimate=affine2d(eye(3));
                flag=0;
            end
            
        
        else  % for multimodal images, imregcorr is more robust 
            if sum(tempx(:))>0 && sum(tempy(:))>0  % ignore blank tiles at the
                tformEstimate=imregtform(tempy,tempx,finalregtype,optimizer,metric,'PyramidLevels',args.pyramid);
                if strcmpi(finalregtype,'affine') 
                    if abs(tformEstimate.T(3,1))>=tError ||abs(tformEstimate.T(3,2))>=tError
                        fprintf('WARNING: Potential failure, displacement = %d,%d pixels. Retrying with rigid.\n',...
                            round(tformEstimate.T(3,1)),round(tformEstimate.T(3,2)));

                        tformEstimate=imregtform(tempy,tempx,'rigid',optimizer,metric,'PyramidLevels',args.pyramid);
                        if abs(tformEstimate.T(3,1))>=tError ||abs(tformEstimate.T(3,2))>=tError
                            fprintf('WARNING: Second potential failure, displacement = %d,%d pixels. Retrying with translation.\n',...
                                round(tformEstimate.T(3,1)),round(tformEstimate.T(3,2)));

                            tformEstimate=imregtform(tempy,tempx,'translation',optimizer,metric,'PyramidLevels',args.pyramid);
                            if abs(tformEstimate.T(3,1))>=tError ||abs(tformEstimate.T(3,2))>=tError
                                fprintf('WARNING: Third potential failure, displacement = %d,%d pixels. Using no registration.\n',...
                                    round(tformEstimate.T(3,1)),round(tformEstimate.T(3,2)));
                                tformEstimate=affine2d(eye(3));

                            end
                        end
                    end
                    
                elseif strcmpi(finalregtype,'rigid') 
                    if abs(tformEstimate.T(3,1))>=tError ||abs(tformEstimate.T(3,2))>=tError
                        fprintf('WARNING: Potential failure, displacement = %d,%d pixels. Retrying with tranlation.\n',...
                            round(tformEstimate.T(3,1)),round(tformEstimate.T(3,2)));

                        tformEstimate=imregtform(tempy,tempx,'translation',optimizer,metric,'PyramidLevels',args.pyramid);
                        
                        if abs(tformEstimate.T(3,1))>=tError ||abs(tformEstimate.T(3,2))>=tError
                            fprintf('WARNING: Second potential failure, displacement = %d,%d pixels. Using no registration.\n',...
                                round(tformEstimate.T(3,1)),round(tformEstimate.T(3,2)));
                            tformEstimate=affine2d(eye(3));

                        end
                    end
                elseif strcmpi(finalregtype,'translation') 
                    if abs(tformEstimate.T(3,1))>=tError ||abs(tformEstimate.T(3,2))>=tError
                        fprintf('WARNING: Potential failure, displacement = %d,%d pixels. Using no registration.\n',...
                                round(tformEstimate.T(3,1)),round(tformEstimate.T(3,2)));
                            tformEstimate=affine2d(eye(3));
                    end
                end

                flag=1;
            else
                tformEstimate=affine2d(eye(3));
                flag=0;
            end
        end
        Rfixed = imref2d(size(tempx));
        tempynew=imwarp(tempy,tformEstimate,'OutputView',Rfixed);
        tempyorig=imwarp(tempyorig,tformEstimate,'OutputView',Rfixed);
        %%
        tform(:,:, count)=tformEstimate.T;
        

        regY{r,c}=tempyorig;
        if flag
            if corr2(single(tempx),single(tempy))<corr2(single(tempx),single(tempynew))
                fprintf('Chunk %d x %d : translation = %.2f x %.2f pixels, corr %.4f-->%.4f\n',...
                    r,c,tformEstimate.T(3,1),tformEstimate.T(3,2),corr2(single(tempx),single(tempy)),...
                    corr2(single(tempx),single(tempynew)));
            else
                fprintf('Chunk %d x %d : translation = %.2f x %.2f pixels, corr %.4f<--%.4f\n',...
                    r,c,tformEstimate.T(3,1),tformEstimate.T(3,2),corr2(single(tempx),single(tempy)),...
                    corr2(single(tempx),single(tempynew)));
            end
        end
        count=count+1;
    end
end

% L=prod(chunks);
% T1=zeros(L,3);
% for i=1:L
%     T1(i,1:2)=tform(3,1:2,i);    
% end
% [X,Y] = meshgrid(1:chunks(1),1:chunks(2));
% figure;quiver(Y,X,reshape(T1(:,1),chunks),reshape(T1(:,2),chunks))
% axis([0 chunks(1)+1 0 chunks(2)+1])


fprintf('Merging transformed chunks.\n');
MergeY = [];
for r=progress(1:chunks(1))
    MergeX = [];
    for c=1:chunks(2)                        
        MergeX = blend(MergeX, single(regY{r,c}),2, 2*padsize, 3);
    end
    MergeY = blend(MergeY, MergeX, 1, 2*padsize, 3);    
end
MergeY=MergeY(padsize+1:end-padsize,padsize+1:end-padsize);
MergeY=uint16(MergeY);
finaldim=size(MergeY);
outvol=zeros(dim,'uint16');
outvol(1:finaldim(1),1:finaldim(2))=MergeY;

        
if ischar(img2)
    [~,s,~]=fileparts(img2);
elseif isnumeric(img2)
    s='Registration';
end


img2_reg=outvol(1:dim1(1),1:dim1(2));
opts.tform=tform;
opts.tformInit=tformInit;
opts.chunks=chunks;
opts.dsfactor=dsfactor;
opts.padsize=padsize;

if write
    s1=[s '_tranformation_matrix.mat'];
    s1=fullfile(outputdir,s1);
    fprintf('Writing %s\n',s1);
    save(s1,'tform','tformInit','chunks','dsfactor','padsize');
    
    s1=[s '_finalreg.tif'];
    s1=fullfile(outputdir,s1);
    fprintf('Writing %s\n',s1);
    % if image dimension is > 2GB, use saveastiff
    if dim1(1)*dim1(2)>=2*(1024^3)
        options.color     = false;
        options.compress  = 'off';
        options.message   = true;
        options.append    = false;
        options.overwrite = true;
        options.big       = true;
        saveastiff(img2_reg,s1,options);
    else
        imwrite(img2_reg,s1,'Compression','deflate');
    end
end
    
toc   


end


function dataC=blend(dataA, dataB, dim, overlap, option)
% dataC=blend(dataA, dataB, dim, overlap, option)
% 
% dim = 1 --> row-wise blend
% dim = 2 --> column-wise blend
% use option = 3 for linear blend
if isempty(dataA)
    dataC = dataB;
elseif dim == 2 % column-wise blend, horizontal
       ny = size(dataA,1);
       temp_start = dataA(:, end-overlap+1: end);
       temp_end = dataB(:, 1:overlap);
    if option == 1   % simple average   
       temp = (temp_start + temp_end)/2;       
    elseif option == 2   % simple replacement 
       temp = [temp_start(:,1:overlap/2), temp_start(:, overlap/2+1:end)];
    elseif option == 3  % linear blend
       ramp = linspace(0, 1, overlap);
       ramp_start = repmat(ramp, [ny, 1]);
       temp = temp_start.*(1- ramp_start) + temp_end.* ramp_start;    
    end
       
     dataC = [dataA(:, 1: end-overlap), temp, dataB(:, overlap+1:end)];
     
  elseif dim == 1  % row-wise blend
   
       nx = size(dataA,2);
       temp_start = dataA(end-overlap+1:end,:);
       temp_end = dataB(1:overlap,:);
    if option == 1   % simple average   
       temp = (temp_start + temp_end)/2;       
    elseif option == 2   % simple replacement 
       temp = [temp_start(1:overlap/2,:); temp_start(overlap/2+1:end,:)];
    elseif option == 3  % linear blend
       ramp = linspace(0, 1, overlap)';
       ramp_start = repmat(ramp, [1, nx]);     
       temp = temp_start.*(1- ramp_start) + temp_end.* ramp_start;    
    end
       
     dataC = [dataA(1: end-overlap,:); temp; dataB(overlap+1:end,:)];
     
end
end





