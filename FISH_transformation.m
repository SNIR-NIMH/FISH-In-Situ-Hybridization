function FISH_transformation(img1,img2,transformation_mat,outputname)
% 
% FISH_transformation(FIXED, MOVING, TFORM, OUTPUT)
% 
% FIXED      Fixed image used for registration. It is also called "reference"
%            image. The output image will have same dimension as this one. It is
%            the first argument of the FIST_registration code.
% MOVING     Moving image (input) to be transformed to the fixed image
% TFORM      Transformation matrix obtained from the FISH_registration code
% OUTPUT     Output (.tif) image, transformed version of moving image


warning off;
[~,~,ext]=fileparts(outputname);
if ~strcmpi(ext,'.tif')  && ~strcmpi(ext,'.tiff')
    outputname=[outputname '.tif'];
end
if isfile(outputname)
    fprintf('WARNING: Output file (%s) exists. I will overwrite.\n',outputname);
    
end
if ischar(img1)
    fprintf('Reading %s\n',img1);
    x=imread(img1);
elseif isnumeric(img1)
    x=single(img1);
end

if ischar(img2)
    fprintf('Reading %s\n',img2);
    y=imread(img2);
elseif isnumeric(img2)
    y=single(img2);
end
dim1=size(x);
dim2=size(y);

if ischar(transformation_mat)
    mat=load(transformation_mat);
elseif isstruct(transformation_mat)
    mat=transformation_mat;
end

dim=[max(dim1(1),dim2(1)) max(dim1(2),dim2(2))];
x1=zeros(dim,'uint16');
y1=zeros(dim,'uint16');
% Padding images to have same size
x1(1:dim1(1),1:dim1(2))=x;
y1(1:dim2(1),1:dim2(2))=y;
x=x1;y=y1;clear x1 y1;

Rfixed = imref2d(size(x));
y = imwarp(y,mat.tformInit,'OutputView',Rfixed);
chunks=mat.chunks;


padsize=mat.padsize;  % load from mat because padsize can vary with input overlap
% padsize=max(round(0.05*(dim./chunks)));  % a 5% overlap between tiles
pdim=floor(dim./chunks)+2*padsize;  % dimension of the overlapping tiles to be registered

fprintf('Padding size = %d pixels.\n',padsize);
fprintf('Registration chunk size = %d x %d pixels.\n',pdim(1),pdim(2));
d=floor(dim./chunks);
regY=cell(chunks(1),chunks(2));
fprintf('Transforming.\n');
count=1;
for r=1:chunks(1)
    for c=1:chunks(2)
        tempx=zeros(pdim,'uint16');
        tempy=zeros(pdim,'uint16');
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
        tempx(I11:I11+adim(1)-1,J11:J11+adim(2)-1)=a;
        
        a=y(I1:I2,J1:J2);
        adim=size(a);
        tempy(I11:I11+adim(1)-1,J11:J11+adim(2)-1)=a;
        
        tf=affine2d(mat.tform(:,:,count));
        Rfixed = imref2d(size(tempx));
        tempynew=imwarp(tempy,tf,'OutputView',Rfixed);
     
        regY{r,c}=tempynew;
        
        count=count+1;
    end
end


fprintf('Stitching.\n');
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
outvol=outvol(1:dim1(1),1:dim1(2));

       

fprintf('Writing %s\n',outputname);
if dim1(1)*dim1(2)>=2*(1024^3)
    options.color     = false;
    options.compress  = 'adobe';
    options.message   = true;
    options.append    = false;
    options.overwrite = true;
    options.big       = true;
    saveastiff(outvol,outputname,options);
else
    imwrite(outvol,outputname,'Compression','deflate');
end
% imwrite(outvol,outputname);


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

