classdef TruncateToCylinderProcessor < mStir.ImageDataProcessor
%     Class for the image filter that zeroes the image outside the cylinder
%     of the same xy-diameter and z-size as those of the image.
    methods
        function self = TruncateToCylinderProcessor()
            self.name = 'TruncateToCylindricalFOVImageProcessor';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_strictly_less_than_radius(self, flag)
%         Specifies whether the area not affected by filtering is strictly
%         inside the cylinder (flag = True) or not (flag = False).
            if flag
                str = 'true';
            else
                str = 'false';
            end
            mStir.setParameter(self.handle,...
                'TruncateToCylindricalFOVImageProcessor',...
                'strictly_less_than_radius', str, 'c')
        end
        function flag = get_strictly_less_than_radius(self)
%         Returns the answer to the question: Is the area not affected by 
%         filtering strictly inside the cylinder?
            flag = mStir.parameter(self.handle,...
                'TruncateToCylindricalFOVImageProcessor',...
                'strictly_less_than_radius', 'i');
        end
    end
end