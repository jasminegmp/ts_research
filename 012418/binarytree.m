%% Classes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef binarytree
% write a description of the class here.
   properties
   % define the properties of the class here, (like fields of a struct)
       left;
       right;
       ts_seg;
       index;
       level;
   end
   methods
   % methods, including the constructor are defined in this block
       function obj = binarytree(ts_seg,index,left,right, level)
       % class constructor
           if(nargin > 0)
             obj.ts_seg  = ts_seg;
             obj.index = index;
             obj.left  = left;
             obj.right = right;
             obj.level = level;
           end
       end
       
%        function obj = buildtree(obj, root, lc_root, rc_root)
%            root.left = lc_root;
%            root.right = rc_root;
% 
%            obj.day = obj.day + numdays;
%        end
   end
end