classdef AcquisitionModel < handle
%     Class for a PET acquisition model that relates an image x to the
%     acquisition data y as
%     (F)    y = [1/n](G x + [a]) + [b]
%     where:
%     G is the geometric (ray tracing) projector from the image voxels
%     to the scanner's pairs of detectors (bins);
%     a and b are otional additive and background terms representing
%     the effects of noise and scattering; assumed to be 0 if not present;
%     n is an optional bin normalization term representing the inverse of
%     detector (bin) efficiencies; assumed to be 1 if not present.
%     The computation of y for a given x by the above formula (F) is
%     referred to as forward projection, and the computation of
%     (B)    z = G' m y
%     where G' is the transpose of G, is referred to as backward projection.
    properties
        name
        handle
    end
    methods
        function self = AcquisitionModel()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_additive_term(self, at)
%         Sets the additive term a in (F);
%         at:  an AcquisitionData object containing a.
            mStir.setParameter(self.handle, 'AcquisitionModel', ...
                'additive_term', at, 'h');
        end
        function set_background_term(self, bt)
%         Sets the background term b in (F);
%         bt:  an AcquisitionData object containing b.
            mStir.setParameter(self.handle, 'AcquisitionModel', ...
                'additive_term', bt, 'h');
        end
        function set_normalisation(self, bin_eff)
%         Sets the normalisation n in (F);
%         bin_eff:  an AcquisitionData object containing bin efficiencies
%                   (the inverse of n).
            mStir.setParameter(self.handle, 'AcquisitionModel', ...
                'additive_term', bin_eff, 'h');
        end
        function set_up(self, acq_templ, img_templ)
%         Prepares this object for performing forward and backward
%         projections;
%         acq_templ:  an AcquisitionData object used as a template for
%                     creating an AcquisitionData object to store forward
%                     projection;
%         img_templ:  an ImageData object used as a template for creating an
%                     ImageData object to store backward projection.
            h = calllib...
                ('mstir', 'mSTIR_setupAcquisitionModel',...
                self.handle, acq_templ.handle, img_templ.handle);
            mUtil.checkExecutionStatus([self.name ':set_up'], h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function ad = forward(self, image, filename)
%         Returns the forward projection of x given by (F);
%         image   :  an ImageData object containing x;
%         filename:  an optional name of the file to store projection data;
%                    if not present, projection data is stored in memory
%                    (not recommended as it can be huge).
            if nargin < 3
                filename = '';
            end
            ad = mStir.AcquisitionData();
            ad.handle = calllib('mstir', 'mSTIR_acquisitionModelFwd',...
                self.handle, image.handle, filename);
            mUtil.checkExecutionStatus([self.name ':forward'], ad.handle)
        end
        function image = backward(self, ad)
%         Returns the backward projection of y giben by (B);
%         ad:  an AcquisitionData object containing y.
            image = mStir.ImageData();
            image.handle = calllib('mstir', 'mSTIR_acquisitionModelBwd',...
                self.handle, ad.handle);
            mUtil.checkExecutionStatus...
                ([self.name ':backward'], image.handle)
        end
    end
end