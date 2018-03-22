/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

This file collects all testjobs relevant for the tissue
parameter structs and the xml parser filling them.

================================================ */



#include <string>
#include "tissueparameters.h"
#include <stdio.h>
#include <iostream>


#define XML_TEST_PATH "Testing/TestData/test_TissueParameters_XML.xml" 





bool test_allocate_MRTissueParameter_successful(void);
bool test_allocate_PETTissueParameter_successful(void);
bool test_allocate_TissueParameter_successful(void);
bool test_read_TissueParameter_label_from_xml( std::string const xml_filepath);
