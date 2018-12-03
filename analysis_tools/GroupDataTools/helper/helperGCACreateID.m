function [ID ] = helperGCACreateID(path)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if isa(path,'MovieData'); 
 add = 1;  
 pathC = path.outputDirectory_;
else 
    pathC = path;
    add = 0; 
end 
  [~,date] = upDirectory(pathC,2 +add,1);
            [~,num] = upDirectory(pathC,1+add,1);
            [~,name] = upDirectory(pathC,3+add,1);
         ID= [name '_' date '_' num];

end

