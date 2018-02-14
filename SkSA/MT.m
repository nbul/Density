  rng('shuffle');
    X = randi(shortside, MTnumber,1);
    Y = randi(longside, MTnumber,1);
    %% Generate angles with a given distribution
    angles = normrnd(0,distribution, [MTnumber 1]);
    for i=1:MTnumber
        if angles(i, 1) <= -90
            angles(i, 1) = angles(i, 1) + 180;
        elseif angles(i, 1) > 90
            angles(i, 1) = angles(i, 1) - 180;
        end
    end
    
    if bundling > 0
        for i=1:ceil(MTnumber*bundling/100)
            X(MTnumber-ceil(MTnumber*bundling/100)+i)=X(i);
            Y(MTnumber-ceil(MTnumber*bundling/100)+i)=Y(i);
            angles(MTnumber-ceil(MTnumber*bundling/100)+i,1)=angles(i,1);
        end
    end
    
    %% Line parameters and start/end points
    a = 1./tand(angles);
    b = Y - a.*X;
    intersect=zeros(MTnumber,4);
    l=zeros(MTnumber,1);
    for i=1:MTnumber
        l(i)=0;
        X_temp = 1;
        Y_temp = ceil(b(i)+a(i));
        if Y_temp>0 && Y_temp<=longside
            l(i)=1;
            intersect(i,1) = X_temp;
            intersect(i,2) = Y_temp;
        end
        
        X_temp = ceil(1 - b(i)/a(i));
        Y_temp = 1;
        if X_temp>0 && X_temp<=shortside
            if l(i) == 0
                l(i) = 1;
                intersect(i,1) = X_temp;
                intersect(i,2) = Y_temp;
            else
                l(i) = l(i)+1;
                intersect(i,3) = X_temp;
                intersect(i,4) = Y_temp;
            end
        end
        
        X_temp = ceil((longside - b(i))/a(i));
        Y_temp = longside;
        if X_temp>0 && X_temp<=shortside
            if l(i) == 0
                l(i) = 1;
                intersect(i,1) = X_temp;
                intersect(i,2) = Y_temp;
            else
                l(i) = l(i)+1;
                intersect(i,3) = X_temp;
                intersect(i,4) = Y_temp;
            end
        end
        
        X_temp = shortside;
        Y_temp = a(i)*shortside +b(i);
        if Y_temp>0 && Y_temp<=longside
            if l(i) == 0
                l(i) = 1;
                intersect(i,1) = X_temp;
                intersect(i,2) = Y_temp;
            else
                l(i) = l(i)+1;
                intersect(i,3) = X_temp;
                intersect(i,4) = Y_temp;
            end
        end
    end
    
    %% Draw lines
    image = zeros(longside,shortside);
    image_MT_gray = image;
    for i = 1:MTnumber
        image_MT = insertShape(image,'line',intersect(i,:), 'linewidth', 3, 'Color', [I I I]);
        image_MT_gray = image_MT_gray + image_MT(:,:,1);
    end
    
    image_MT_gray = image_MT_gray + 10;
    image_MT_gray(image_MT_gray>255) = 255;
    image_MT_gray = imgaussfilt(image_MT_gray,1);