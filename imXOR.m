function opImg = imXOR(ipImg)
    
    pdImg = zeros(size(ipImg)+[2 2]); 
    pdImg(2:end-1,2:end-1) = ipImg;
    
    vals1 = xor(ipImg,pdImg(1:end-2,1:end-2))*(2^7);
    vals2 = xor(ipImg,pdImg(1:end-2,2:end-1))*(2^6);
    vals3 = xor(ipImg,pdImg(1:end-2,3:end))*(2^5);
    vals4 = xor(ipImg,pdImg(2:end-1,3:end))*(2^4);
    vals5 = xor(ipImg,pdImg(3:end,3:end))*(2^3);
    vals6 = xor(ipImg,pdImg(3:end,2:end-1))*(2^2);
    vals7 = xor(ipImg,pdImg(3:end,1:end-2))*(2^1);
    vals8 = xor(ipImg,pdImg(2:end-1,1:end-2))*(2^0);

    opImg = vals1 + vals2 + vals3 + vals4 + vals5 + vals6 + vals7 + vals8; 

end
