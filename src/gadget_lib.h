#ifndef GADGETS_LIBRARY
#define GADGETS_LIBRARY

class aGadget {
public:
	//	virtual ~aGadget() {}
	virtual std::string xml() const = 0;
};

class IsmrmrdAcqMsgReader : public aGadget {
public:
	virtual std::string xml() const {
		std::string xml_script("<reader>\n");
		xml_script += " <slot>1008</slot>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += 
			" <classname>GadgetIsmrmrdAcquisitionMessageReader</classname>\n";
		xml_script += "</reader>\n";
		return xml_script;
	}
};

class IsmrmrdImgMsgReader : public aGadget {
public:
	virtual std::string xml() const {
		std::string xml_script("<reader>\n");
		xml_script += " <slot>1022</slot>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script +=
			" <classname>MRIImageReader</classname>\n";
		xml_script += "</reader>\n";
		return xml_script;
	}
};

class IsmrmrdImgMsgWriter : public aGadget {
public:
	virtual std::string xml() const {
		std::string xml_script("<writer>\n");
		xml_script += " <slot>1022</slot>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>MRIImageWriter</classname>\n";
		xml_script += "</writer>\n";
		return xml_script;
	}
};

class RemoveOversamplingGadget : public aGadget {
public:
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
	        xml_script += " <name>RemoveROOversampling</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>RemoveROOversamplingGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class AcqAccTrigGadget : public aGadget {
public:
	AcqAccTrigGadget() : trigger_dimension_("repetition"), sorting_dimension_("slice") {}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
	        xml_script += " <name>AccTrig</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>AcquisitionAccumulateTriggerGadget</classname>\n";
	        xml_script += " <property>\n";
	        xml_script += "  <name>trigger_dimension</name>\n";
	        xml_script += "  <value>" + trigger_dimension_ + "</value>\n";
	        xml_script += " </property>\n";
	        xml_script += " <property>\n";
	        xml_script += "  <name>sorting_dimension</name>\n";
	        xml_script += "  <value>" + sorting_dimension_ + "</value>\n";
	        xml_script += " </property>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
private:
	std::string trigger_dimension_;
	std::string sorting_dimension_;
};

class BucketToBuffGadget : public aGadget {
public:
	BucketToBuffGadget() : n_dimension_(""), s_dimension_(""), split_slices_("true") {}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
	        xml_script += " <name>Buff</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>BucketToBufferGadget</classname>\n";
	        xml_script += " <property>\n";
	        xml_script += "  <name>N_dimension</name>\n";
	        xml_script += "  <value>" + n_dimension_ + "</value>\n";
	        xml_script += " </property>\n";
	        xml_script += " <property>\n";
	        xml_script += "  <name>S_dimension</name>\n";
	        xml_script += "  <value>" + s_dimension_ + "</value>\n";
	        xml_script += " </property>\n";
	        xml_script += " <property>\n";
	        xml_script += "  <name>split_slices</name>\n";
	        xml_script += "  <value>" + split_slices_ + "</value>\n";
	        xml_script += " </property>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
private:
	std::string n_dimension_;
	std::string s_dimension_;
	std::string split_slices_;
};

class SimpleReconstructionGadget : public aGadget {
public:
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
	        xml_script += " <name>SimpleRecon</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>SimpleReconGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class ImgArrSplitGadget : public aGadget {
public:
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
	        xml_script += " <name>ImageArraySplit</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>ImageArraySplitGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class ExtGadget : public aGadget {
public:
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
	        xml_script += " <name>Extract</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>ExtractGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class ImgFinishGadget : public aGadget {
public:
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
	        xml_script += " <name>ImageFinish</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>ImageFinishGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

#endif
