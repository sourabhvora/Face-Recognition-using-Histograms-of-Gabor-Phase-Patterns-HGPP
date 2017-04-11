function [ HGPP ] = compute_HGPP( filtered_image, r, c, num_pixels, num_scales, num_orientations, sub_region_size, num_bins )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%Create Daugman's encoding
P_real = zeros(size(filtered_image));
P_imag = zeros(size(filtered_image));

P_real(real(filtered_image)<=0) = 1;
P_imag(imag(filtered_image)<=0) = 1;


% Compute Global Gabor Phase Patterns
GGPP_real = zeros(num_scales,num_pixels);
GGPP_imag = zeros(num_scales,num_pixels);

%{
tic;
for i=1:num_scales
    for pix=1:num_pixels
        s_real = '';
        s_imag = '';
        for j=0:num_orientations-1
            col_idx = (i-1)*num_pixels*num_orientations + j*num_pixels + pix;
            s_real = [s_real num2str(P_real(col_idx))];
            s_imag = [s_imag num2str(P_imag(col_idx))];
        end
        GGPP_real(i,pix) = bin2dec(s_real);     % i starts from 0. Hence i+1.
        GGPP_imag(i,pix) = bin2dec(s_imag);     % i starts from 0. Hence i+1.
    end
end
toc;
%}

for i=1:num_scales
    s_real=zeros(num_pixels,1);
    s_imag=zeros(num_pixels,1);
    
    col_idx = (i-1)*num_pixels*num_orientations + 0*num_pixels;
    img_real = P_real(col_idx+1:col_idx+num_pixels);
    img_imag = P_imag(col_idx+1:col_idx+num_pixels);
    s_real = s_real + (2^7)*img_real;
    s_imag = s_imag + (2^7)*img_imag;
    
    
    col_idx = (i-1)*num_pixels*num_orientations + 1*num_pixels;
    img_real = P_real(col_idx+1:col_idx+num_pixels);
    img_imag = P_imag(col_idx+1:col_idx+num_pixels);
    s_real = s_real + (2^6)*img_real;
    s_imag = s_imag + (2^6)*img_imag;
    
    
    col_idx = (i-1)*num_pixels*num_orientations + 2*num_pixels;
    img_real = P_real(col_idx+1:col_idx+num_pixels);
    img_imag = P_imag(col_idx+1:col_idx+num_pixels);
    s_real = s_real + (2^5)*img_real;
    s_imag = s_imag + (2^5)*img_imag;
    
    
    col_idx = (i-1)*num_pixels*num_orientations + 3*num_pixels;
    img_real = P_real(col_idx+1:col_idx+num_pixels);
    img_imag = P_imag(col_idx+1:col_idx+num_pixels);
    s_real = s_real + (2^4)*img_real;
    s_imag = s_imag + (2^4)*img_imag;
    
    col_idx = (i-1)*num_pixels*num_orientations + 4*num_pixels;
    img_real = P_real(col_idx+1:col_idx+num_pixels);
    img_imag = P_imag(col_idx+1:col_idx+num_pixels);
    s_real = s_real + (2^3)*img_real;
    s_imag = s_imag + (2^3)*img_imag;
    
    
    col_idx = (i-1)*num_pixels*num_orientations + 5*num_pixels;
    img_real = P_real(col_idx+1:col_idx+num_pixels);
    img_imag = P_imag(col_idx+1:col_idx+num_pixels);
    s_real = s_real + (2^2)*img_real;
    s_imag = s_imag + (2^2)*img_imag;
    
    
    col_idx = (i-1)*num_pixels*num_orientations + 6*num_pixels;
    img_real = P_real(col_idx+1:col_idx+num_pixels);
    img_imag = P_imag(col_idx+1:col_idx+num_pixels);
    s_real = s_real + (2^1)*img_real;
    s_imag = s_imag + (2^1)*img_imag;
    
    
    col_idx = (i-1)*num_pixels*num_orientations + 7*num_pixels;
    img_real = P_real(col_idx+1:col_idx+num_pixels);
    img_imag = P_imag(col_idx+1:col_idx+num_pixels);
    s_real = s_real + (2^0)*img_real;
    s_imag = s_imag + (2^0)*img_imag;
    
    GGPP_real(i,:) = s_real';
    GGPP_imag(i,:) = s_imag';
end


% Compute Local Gabor Phase Patterns
LGPP_real = cell(num_scales, num_orientations);
LGPP_imag = cell(num_scales, num_orientations);

for i=1:num_scales
    for j=1:num_orientations
        col_idx = (i-1)*num_pixels*num_orientations + (j-1)*num_pixels;
        
        img_uv_real = P_real(col_idx+1:col_idx+num_pixels);
        img_uv_real = reshape(img_uv_real,r,c);
        
        img_uv_real = imXOR(img_uv_real);
        
        img_uv_imag = P_imag(col_idx+1:col_idx+num_pixels);
        img_uv_imag = reshape(img_uv_imag,r,c);
        
        img_uv_imag = imXOR(img_uv_imag);
        %{
        for k=2:r-1
            for l=2:c-1
                s = [ num2str(xor(img_uv_real(k,l),img_uv_real(k-1,l-1)))  num2str(xor(img_uv_real(k,l),img_uv_real(k-1,l)))  num2str(xor(img_uv_real(k,l),img_uv_real(k-1,l+1)))  num2str(xor(img_uv_real(k,l),img_uv_real(k,l+1)))  num2str(xor(img_uv_real(k,l),img_uv_real(k+1,l+1)))  num2str(xor(img_uv_real(k,l),img_uv_real(k+1,l)))  num2str(xor(img_uv_real(k,l),img_uv_real(k+1,l-1)))  num2str(xor(img_uv_real(k,l),img_uv_real(k,l-1)))  ];
                img_uv_real(k,l) = bin2dec(s);
                
                s = [ num2str(xor(img_uv_imag(k,l),img_uv_imag(k-1,l-1)))  num2str(xor(img_uv_imag(k,l),img_uv_imag(k-1,l)))  num2str(xor(img_uv_imag(k,l),img_uv_imag(k-1,l+1)))  num2str(xor(img_uv_imag(k,l),img_uv_imag(k,l+1)))  num2str(xor(img_uv_imag(k,l),img_uv_imag(k+1,l+1)))  num2str(xor(img_uv_imag(k,l),img_uv_imag(k+1,l)))  num2str(xor(img_uv_imag(k,l),img_uv_imag(k+1,l-1)))  num2str(xor(img_uv_imag(k,l),img_uv_imag(k,l-1)))  ];
                img_uv_imag(k,l) = bin2dec(s);
            end
        end
        %}
        LGPP_real{i,j} = img_uv_real(:);
        LGPP_imag{i,j} = img_uv_imag(:);
    end
end

% Compute Histogram of Gabor Phase Patterns

HGPP = [];

% 1. H_GGPP_real

for i=1:num_scales
    im = reshape(GGPP_real(i,:),r,c);
    for j=1:sub_region_size(1):r
        for k=1:sub_region_size(1):c
            x = im(j:j+sub_region_size(1)-1,k:k+sub_region_size(2)-1);
            [counts, ~] = imhist(x, num_bins);
            HGPP = [HGPP;counts];
        end
    end
end

% 2. H_GGPP_imag

for i=1:num_scales
    im = reshape(GGPP_imag(i,:),r,c);
    for j=1:sub_region_size(1):r
        for k=1:sub_region_size(1):c
            x = im(j:j+sub_region_size(1)-1,k:k+sub_region_size(2)-1);
            [counts, ~] = imhist(x, num_bins);
            HGPP = [HGPP;counts];
        end
    end
end

% 3. H_LGPP_real

for h=1:num_scales
    for i=1:num_orientations
        im = reshape(LGPP_real{h,i},r,c);
        for j=1:sub_region_size(1):r
            for k=1:sub_region_size(1):c
                x = im(j:j+sub_region_size(1)-1,k:k+sub_region_size(2)-1);
                [counts, ~] = imhist(x, num_bins);
                HGPP = [HGPP;counts];
            end
        end
    end
end

% 4. H_LGPP_imag

for h=1:num_scales
    for i=1:num_orientations
        im = reshape(LGPP_imag{h,i},r,c);
        for j=1:sub_region_size(1):r
            for k=1:sub_region_size(1):c
                x = im(j:j+sub_region_size(1)-1,k:k+sub_region_size(2)-1);
                [counts, ~] = imhist(x, num_bins);
                HGPP = [HGPP;counts];
            end
        end
    end
end




end

