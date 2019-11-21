function ret = fn_read_file(path, file_name, num_col)
    fname = sprintf('%s/%s.txt', path, file_name);
    fileID = fopen(fname,'r');
    
    switch num_col
        case 1
            ret = fscanf(fileID, '%f', [1 inf]);
        case 2
            ret = fscanf(fileID, '%f %f', [2 inf]);
        case 3
            ret = fscanf(fileID, '%f %f %f', [3 inf]);
        case 4
            ret = fscanf(fileID, '%f %f %f %f', [4 inf]);
        case 5
            ret = fscanf(fileID, '%f %f %f %f %f', [5 inf]);
        case 6
            ret = fscanf(fileID, '%f %f %f %f %f %f', [6 inf]);
        case 8
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f', [8 inf]);
        case 9
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f', [9, inf]);
        case 10
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f', [10, inf]);
        case 11
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f', [11, inf]);
        case 12
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f', [12 inf]);
        case 13
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f', [13 inf]);    
        case 15
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [15 inf]);    
        case 20
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [20 inf]);
        case 21
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [21 inf]);
        case 22
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [22 inf]);   
        case 23
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [23 inf]);   
        case 25
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [25 inf]);   
        case 26
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [26 inf]);   
        case 27
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [27 inf]);   
        case 28
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [28 inf]);   
        case 29
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [29 inf]);
        case 33
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [33 inf]);
        case 34
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [34 inf]);           
        case 35
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [35 inf]);           
        case 39
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [39 inf]);           
        case 45
            ret = fscanf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [45 inf]);           
    end
    fclose(fileID);
end