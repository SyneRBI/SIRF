#ifndef GADGETS_LIBRARY
#define GADGETS_LIBRARY

#include <boost/algorithm/string.hpp>

#include "data_handle.h"

class aGadget {
public:
	//	virtual ~aGadget() {}
	std::string name() const {
		return this->name_;
	}
	virtual void set_property(const char* prop, const char* value) = 0;
	virtual std::string xml() const = 0;
protected:
	std::string name_;
};

class IsmrmrdAcqMsgReader : public aGadget {
public:
	IsmrmrdAcqMsgReader() {
		name_ = "GadgetIsmrmrdAcquisitionMessageReader";
	}
	virtual void set_property(const char* prop, const char* value) {}
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

class IsmrmrdAcqMsgWriter : public aGadget {
public:
	IsmrmrdAcqMsgWriter() {
		name_ = "GadgetIsmrmrdAcquisitionMessageWriter";
	}
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string xml() const {
		std::string xml_script("<writer>\n");
		xml_script += " <slot>1008</slot>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script +=
			" <classname>GadgetIsmrmrdAcquisitionMessageWriter</classname>\n";
		xml_script += "</writer>\n";
		return xml_script;
	}
};

class IsmrmrdImgMsgReader : public aGadget {
public:
	IsmrmrdImgMsgReader() {
		name_ = "MRIImageReader";
	}
	virtual void set_property(const char* prop, const char* value) {}
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
	IsmrmrdImgMsgWriter() {
		name_ = "MRIImageWriter";
	}
	virtual void set_property(const char* prop, const char* value) {}
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
	RemoveOversamplingGadget() {
		name_ = "RemoveROOversamplingGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
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
	AcqAccTrigGadget() : trigger_dimension_("repetition"), sorting_dimension_("slice") {
		name_ = "AcquisitionAccumulateTriggerGadget";
	}
	virtual void set_property(const char* prop, const char* value) 
	{
		if (boost::iequals(prop, "trigger_dimension"))
			trigger_dimension_ = value;
		else if (boost::iequals(prop, "sorting_dimension"))
			sorting_dimension_ = value;
		else
			THROW("unknown gadget parameter");
	}
	virtual std::string xml() const 
	{
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
	BucketToBuffGadget() : n_dimension_(""), s_dimension_(""), split_slices_("true") {
		name_ = "BucketToBufferGadget";
	}
	virtual void set_property(const char* prop, const char* value) 
	{
		if (boost::iequals(prop, "n_dimension"))
			n_dimension_ = value;
		else if (boost::iequals(prop, "s_dimension"))
			s_dimension_ = value;
		else if (boost::iequals(prop, "split_slices"))
			split_slices_ = value;
		else
			THROW("unknown gadget parameter");
	}
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
	SimpleReconstructionGadget() {
		name_ = "SimpleReconGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
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
	ImgArrSplitGadget() {
		name_ = "ImageArraySplitGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
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
	ExtGadget() {
		name_ = "ExtractGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
    xml_script += " <name>Extract</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>ExtractGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class ComplexToFloatGadget : public aGadget {
public:
	ComplexToFloatGadget() {
		name_ = "ComplexToFloatGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
		xml_script += " <name>ComplexToFloatAttrib</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>ComplexToFloatGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class FloatToShortGadget : public aGadget {
public:
	FloatToShortGadget()
	{
		name_ = "FloatToShortGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string xml() const 
	{
		std::string xml_script("<gadget>\n");
		xml_script += " <name>FloatToShortAttrib</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>FloatToUShortGadget</classname>\n";
		xml_script += "<property><name>max_intensity</name><value>32767</value></property>";
		xml_script += "<property><name>min_intensity</name><value>0</value></property>";
		xml_script += "<property><name>intensity_offset</name><value>0</value></property>";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class ImgFinishGadget : public aGadget {
public:
	ImgFinishGadget() {
		name_ = "ImageFinishGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
		xml_script += " <name>ImageFinish</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>ImageFinishGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class AcqFinishGadget : public aGadget {
public:
	AcqFinishGadget() {
		name_ = "AcquisitionFinishGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
		xml_script += " <name>AcquisitionFinish</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>AcquisitionFinishGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class SimpleReconstructionGadgetSet : public aGadget {
public:
	SimpleReconstructionGadgetSet() 
	{
		name_ = "SimpleReconGadgetSet";
	}
	virtual void set_property(const char* prop, const char* value) 
	{
		if (boost::iequals(prop, "trigger_dimension"))
			aat_.set_property(prop, value);
		else if (boost::iequals(prop, "sorting_dimension"))
			aat_.set_property(prop, value);
		else if (boost::iequals(prop, "n_dimension"))
			bb_.set_property(prop, value);
		else if (boost::iequals(prop, "s_dimension"))
			bb_.set_property(prop, value);
		else if (boost::iequals(prop, "split_slices"))
			bb_.set_property(prop, value);
		else
			THROW("unknown gadget parameter");
	}
	virtual std::string xml() const 
	{
		std::string xml_script;
		xml_script += aat_.xml();
		xml_script += bb_.xml();
		xml_script += sr_.xml();
		xml_script += ias_.xml();
		return xml_script;
	}
private:
	AcqAccTrigGadget aat_;
	BucketToBuffGadget bb_;
	SimpleReconstructionGadget sr_;
	ImgArrSplitGadget ias_;
	//std::string trigger_dimension_;
	//std::string sorting_dimension_;
	//std::string n_dimension_;
	//std::string s_dimension_;
	//std::string split_slices_;
};

#endif
