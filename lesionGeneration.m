function [img] = lesionGeneration(img, lesioncoef, lesionsize)
ROI = size(img,1);
zSlice = size(img,3);
resultimg = zeros(ROI, ROI, ROI);
resultimg(ceil(ROI/2),ceil(ROI/2),ceil(ROI/2)) = 1;
directions = ['L', 'R', 'F', 'B', 'U', 'D'];
for pixel_num = 1:lesionsize
    direction = directions(randi([1 6]));
        resultimg = growthlesion(resultimg, direction);
    
end
center_slice = floor(ROI/2+1);
resultimg = resultimg(:,:,center_slice-floor(zSlice/2):center_slice+floor(zSlice/2));
img(resultimg == 1) = lesioncoef;
end


function [resultImg] = growthlesion(targetImg, direction)
datasize = size(targetImg,1 );
switch direction
    case 'L'
        trans = max(targetImg, [], 1);
        axisIndex = find(trans == 1);
        targetIndex = datasample(axisIndex, 1);
        targetIndex1 = mod(targetIndex, datasize);
        targetIndex2 = (targetIndex - targetIndex1)/datasize + 1;
        pixelList = targetImg(:,targetIndex1, targetIndex2);
        loopFlag = 0;
        idx = 1;
        while loopFlag == 0
            if pixelList(idx) == 0
                idx = idx + 1;
            else
                loopFlag = 1;
            end
        end
        resultImg = targetImg;
        resultImg(idx - 1, targetIndex1, targetIndex2) = 1;
    case 'R'
        trans = max(targetImg, [], 1);
        axisIndex = find(trans == 1);
        targetIndex = datasample(axisIndex, 1);
        targetIndex1 = mod(targetIndex, datasize);
        targetIndex2 = (targetIndex - targetIndex1)/datasize + 1;
        pixelList = targetImg(:,targetIndex1, targetIndex2);
        loopFlag = 0;
        idx = datasize;
        while loopFlag == 0
            if pixelList(idx) == 0
                idx = idx - 1;
            else
                loopFlag = 1;
            end
        end
        resultImg = targetImg;
        resultImg(idx + 1, targetIndex1, targetIndex2) = 1;
    case 'F'
        trans = max(targetImg, [], 2);
        axisIndex = find(trans == 1);
        targetIndex = datasample(axisIndex, 1);
        targetIndex1 = mod(targetIndex, datasize);
        targetIndex2 = (targetIndex - targetIndex1)/datasize + 1;
        pixelList = targetImg(targetIndex1, :, targetIndex2);
        loopFlag = 0;
        idx = 1;
        while loopFlag == 0
            if pixelList(idx) == 0
                idx = idx + 1;
            else
                loopFlag = 1;
            end
        end
        resultImg = targetImg;
        resultImg( targetIndex1, idx - 1, targetIndex2) = 1;
    case 'B'
        trans = max(targetImg, [], 2);
        axisIndex = find(trans == 1);
        targetIndex = datasample(axisIndex, 1);
        targetIndex1 = mod(targetIndex, datasize);
        targetIndex2 = (targetIndex - targetIndex1)/datasize + 1;
        pixelList = targetImg(targetIndex1, :, targetIndex2);
        loopFlag = 0;
        idx = datasize;
        while loopFlag == 0
            if pixelList(idx) == 0
                idx = idx - 1;
            else
                loopFlag = 1;
            end
        end
        resultImg = targetImg;
        resultImg( targetIndex1, idx + 1, targetIndex2) = 1;
    case 'U'
        trans = max(targetImg, [], 3);
        axisIndex = find(trans == 1);
        targetIndex = datasample(axisIndex, 1);
        targetIndex1 = mod(targetIndex, datasize);
        targetIndex2 = (targetIndex - targetIndex1)/datasize + 1;
        pixelList = targetImg(targetIndex1, targetIndex2, :);
        loopFlag = 0;
        idx = 1;
        while loopFlag == 0
            if pixelList(idx) == 0
                idx = idx + 1;
            else
                loopFlag = 1;
            end
        end
        resultImg = targetImg;
        resultImg( targetIndex1,targetIndex2,  idx - 1) = 1;
    case 'D'
        trans = max(targetImg, [], 3);
        axisIndex = find(trans == 1);
        targetIndex = datasample(axisIndex, 1);
        targetIndex1 = mod(targetIndex, datasize);
        targetIndex2 = (targetIndex - targetIndex1)/datasize + 1;
        pixelList = targetImg(targetIndex1, targetIndex2, :);
        loopFlag = 0;
        idx = datasize;
        while loopFlag == 0
            if pixelList(idx) == 0
                idx = idx - 1;
            else
                loopFlag = 1;
            end
        end
        resultImg = targetImg;
        resultImg( targetIndex1,targetIndex2,  idx + 1) = 1;
end
end
