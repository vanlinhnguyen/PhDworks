function [isLength isUnlimited] = The_dimensions(ncInfo,dim_name)

%initialisation 
isLength=[];
isUnlimited=[];

found = 1;%=0 if dimension not found

%position in ncInfo (ncInfo = nc.getInfo; where nc is a mDataset object)
position = 1;
%check in ncInfo where is dim_name
dim_to_find = '';
while(~strcmp(dim_to_find,dim_name))
    dim_to_find=ncInfo(position:position+length(dim_name)-1);
    position=position+1;%next character
    if position ==(length(ncInfo)-length(dim_name))
        disp('dimension not found')
        found = 0;
    break
    end
end

%position is on the second character of dim_name

if (found == 1)%dimension exist
    position = position + length(dim_name) + 2; %position is now at the begining of its value (+2 for =' ')
    if (ncInfo(position) == 'U')%this dimension is unlimited
        isUnlimited = 'yes';
        position = position + 17;%17 characters for "UNLIMITED;   // (." where . is the first character of the dimension length
        isLength = '';
        number_in_isLength = 1;
        while (ncInfo(position)~=' ')%the dimension length is finish when we get the space character
            isLength(number_in_isLength) = ncInfo(position);%get all the numbers in isLength
            number_in_isLength = number_in_isLength +1;
            position = position + 1;
        end
        isLength = str2num(isLength);%convert to number
    else%not unlimited
        isUnlimited = 'no';
        isLength = '';
        number_in_isLength = 1;
        while (ncInfo(position)~=' ')%the dimension length is finish when we get the space character
            isLength(number_in_isLength) = ncInfo(position);%get all the numbers in isLength
            number_in_isLength = number_in_isLength +1;
            position = position + 1;
        end
        isLength = str2num(isLength);%convert to number  
    end
end