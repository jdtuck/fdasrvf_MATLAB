function resample_spline3()
A = imread('moon.tif');
q = pi/4;
A = A(200:400,100:300);
figure; imshow(A);
resamp = makeresampler('Type','custom','Ndims',Inf, 'ResampleFcn', @imtb_resample_spline3);
resamp2 = makeresampler('cubic','symmetric');
tform = maketform('affine',[cos(q) sin(q) 0; -sin(q) cos(q) 0; 0 0 1]);
B = imtransform(A,tform,resamp);
B2 = imtransform(A,tform,resamp2);
figure; subplot(1,2,1); imshow(B/255); subplot(1,2,2); imshow(B2);
end

function B = imtb_resample_spline3(A,M,TDIMS_A,TDIMS_B,FSIZE_A,FSIZE_B,F,R)
% Spline 3 resampler which supports the image processing toolbox interface 
% resample_fcn(A,M,TDIMS_A,TDIMS_B,FSIZE_A,FSIZE_B,F,R)
    if(nargin<4)
        B = [];
        return;
    end

    A = permute(A, TDIMS_A);

    Asp = Bspline(double(A),3);
    
    Mnat = num2cell(M, [1,2]);
    
    B = Asp(Mnat{:});
    
    B = permute(B, TDIMS_B);
end