classdef MRReconstruction < handle
    properties
        connector_
        images_
    end
    methods
        function self = MRReconstruction()
            self.connector_ = gadgetron.ClientConnector();
            self.images_ = gadgetron.ImagesList();
            self.connector_.register_images_receiver(self.images_);
        end
        function process(self, gc, input_data)
            self.connector_.connect('localhost', '9002')
            self.connector_.config_gadget_chain(gc)
            self.connector_.send_parameters(input_data.header_)
            self.connector_.send_acquisitions(input_data)
            self.connector_.disconnect()
        end
%         function process(self, conn, gc, input_data)
%             self.images_ = gadgetron.ImagesList();
%             conn.register_images_receiver(self.images_);
%             conn.connect('localhost', '9002')
%             conn.config_gadget_chain(gc)
%             conn.send_parameters(input_data.header_)
%             conn.send_acquisitions(input_data)
%             conn.disconnect()
%         end
        function images = get_output(self)
            images = gadgetron.ImagesList(self.images_);
        end
    end
end