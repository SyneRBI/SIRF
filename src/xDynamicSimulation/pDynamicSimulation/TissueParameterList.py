import xmltodict 
import numpy as np


class TissueParameterList:

    def __init__(self):

        self.tissue_parameters = np.array([], dtype=object)

    def print_contents(self):
        for tp in self.tissue_parameters:
            msg = "We have extracted {} for label {}.".format(tp.name, tp.label)
            print(msg)

    def read_xml_to_dict(xml_path):

        with open(xml_path, 'r') as file:
            xml_data = file.read()
        xml_dict = xmltodict.parse(xml_data)
    
        return xml_dict

    def parse_xml(self, xml_path):

        xml_dict = TissueParameterList.read_xml_to_dict(xml_path)
        self.extract_dictionary(xml_dict)
    

    def extract_dictionary(self, xml_dict):

        parameter_list = xml_dict['TissueParameterList']['TissueParameter']

        for param in parameter_list:
            label = int(param['label'])
            name = param['name']

            mr_param = param['MRTissueParameter']

            rho   = float(mr_param['spin_density_percentH2O'])
            t1    = float(mr_param['t1_miliseconds'])
            t2    = float(mr_param['t2_miliseconds'])
            cs    = float(mr_param['cs_ppm'])
            
            pet_param = param['PETTissueParameter']
            mu = float(pet_param['attenuation_1_by_cm'])
            act = float(pet_param['activity_kBq_ml'])

            tp = self.TissueParameter(label, name, rho, t1, t2, cs, mu, act)
            self.tissue_parameters = np.append(self.tissue_parameters, tp)
        
    def get_labels_and_names(self):
        
        label_dict = {}

        for tp in self.tissue_parameters:
            label_dict[tp.name] = tp.label

        return label_dict

    def mr_as_array(self):

        arr = np.empty( shape=(self.tissue_parameters.size, 5))

        for i in range(self.tissue_parameters.size):
            arr[i,:] = self.tissue_parameters[i].mr_as_array()
        
        return arr

    class TissueParameter:

        def __init__(self, label, name, spin_density, T1_ms, T2_ms, cs_ppm, activity_kBq_ml, attenuation_1_by_cm):

            self.label = label
            self.name = name
            self.spin_density = spin_density
            self.T1_ms = T1_ms
            self.T2_ms = T2_ms
            self.cs_ppm = cs_ppm
            self.activity_kBq_ml = activity_kBq_ml
            self.attenuation_1_by_cm = attenuation_1_by_cm
        
        def as_array(self):
            return np.array([self.label, self.spin_density, self.T1_ms, self.T2_ms, self.cs_ppm, self.attenuation_1_by_cm, self.activity_kBq_ml])
            
        def mr_as_array(self):
            arr = self.as_array()
            extracted_idx = [0, 1, 2, 3, 4]
            return arr[extracted_idx]
        
        def pet_as_array(self):
            arr = self.as_array()
            extracted_idx = [0, 5, 6]
            return arr[extracted_idx]
