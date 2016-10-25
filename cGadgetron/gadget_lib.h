#ifndef GADGETS_LIBRARY
#define GADGETS_LIBRARY

#include <map>
#include <boost/algorithm/string.hpp>

class aGadget {
public:
	virtual void set_property(const char* prop, const char* value) = 0;
	virtual std::string value_of(const char* prop) = 0;
	virtual std::string xml() const = 0;
};

class Gadget : public aGadget {
public:
	Gadget(std::string name, std::string dll, std::string cl) :
		gadget_(name), dll_(dll), class_(cl)
	{	}
	virtual void set_property(const char* prop, const char* value)
	{
		par_[prop] = value;
	}
	virtual std::string value_of(const char* prop)
	{
		return par_[prop];
	}
	void add_property(const char* prop, const char* value)
	{
		par_[prop] = value;
	}
	virtual std::string xml() const
	{
		std::string xml_script("<gadget>\n");
		xml_script += " <name>" + gadget_ + "</name>\n";
		xml_script += " <dll>" + dll_ +"</dll>\n";
		xml_script += " <classname>" + class_ + "</classname>\n";
#ifdef MSVC
		std::map<std::string, std::string>::const_iterator it;
#else
		typename std::map<std::string, std::string>::const_iterator it;
#endif
		for (it = par_.begin(); it != par_.end(); ++it) {
			xml_script += " <property>\n";
			xml_script += "  <name>" + it->first + "</name>\n"; \
			xml_script += "  <value>" + it->second + "</value>\n"; \
			xml_script += " </property>\n";
		}
		xml_script += "</gadget>\n";
		return xml_script;
	}
protected:
	std::string gadget_;
	std::string dll_;
	std::string class_;
	std::map<std::string, std::string> par_;
};

class IsmrmrdAcqMsgReader : public aGadget {
public:
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string value_of(const char* prop)
	{
		return std::string("");
	}
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
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string value_of(const char* prop)
	{
		return std::string("");
	}
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
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string value_of(const char* prop)
	{
		return std::string("");
	}
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
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string value_of(const char* prop)
	{
		return std::string("");
	}
	virtual std::string xml() const {
		std::string xml_script("<writer>\n");
		xml_script += " <slot>1022</slot>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>MRIImageWriter</classname>\n";
		xml_script += "</writer>\n";
		return xml_script;
	}
};

class NoiseAdjustGadget : public Gadget {
public:
	NoiseAdjustGadget() :
		Gadget("NoiseAdjust", "gadgetron_mricore", "NoiseAdjustGadget")
	{}
};

class AsymmetricEchoGadget : public Gadget {
public:
	AsymmetricEchoGadget() :
		Gadget("AsymmetricEcho", "gadgetron_mricore", 
		"AsymmetricEchoAdjustROGadget")
	{}
};

class RemoveOversamplingGadget : public Gadget {
public:
	RemoveOversamplingGadget() :
		Gadget("RemoveROOversampling", "gadgetron_mricore",
		"RemoveROOversamplingGadget")
	{}
};

class AcqAccTrigGadget : public Gadget {
public:
	AcqAccTrigGadget() :
		Gadget("AccTrig", "gadgetron_mricore", "AcquisitionAccumulateTriggerGadget")
	{
		add_property("trigger_dimension", "repetition");
		add_property("sorting_dimension", "slice");
	}
};

class BucketToBuffGadget : public Gadget {
public:
	BucketToBuffGadget() :
		Gadget("Buff", "gadgetron_mricore", "BucketToBufferGadget")
	{
		add_property("N_dimension", "");
		add_property("S_dimension", "");
		add_property("split_slices", "true");
		add_property("ignore_segment", "true");
		add_property("verbose", "true");
	}
};

class PrepRefGadget : public Gadget {
public:
	PrepRefGadget() :
		Gadget("PrepRef", "gadgetron_mricore", 
		"GenericReconCartesianReferencePrepGadget")
	{
		add_property("debug_folder", "");
		add_property("perform_timing", "true");
		add_property("verbose", "true");
		add_property("average_all_ref_N", "true");
		add_property("average_all_ref_S", "true");
		add_property("prepare_ref_always", "true");
	}
};

class SimpleReconstructionGadget : public Gadget {
public:
	SimpleReconstructionGadget() :
		Gadget("SimpleRecon", "gadgetron_mricore", "SimpleReconGadget")
	{}
};

class CartesianGrappaGadget : public Gadget {
public:
	CartesianGrappaGadget() :
		Gadget("CartesianGrappa", "gadgetron_mricore",
		"GenericReconCartesianGrappaGadget")
	{
		add_property("image_series", "0");
		add_property("coil_map_algorithm", "Inati");
		add_property("downstream_coil_compression", "true");
		add_property("downstream_coil_compression_thres", "0.01");
		add_property("downstream_coil_compression_num_modesKept", "0");
		add_property("send_out_gfactor", "true");
		add_property("debug_folder", "");
		add_property("perform_timing", "true");
		add_property("verbose", "true");
	}
};

class FOVAdjustmentGadget : public Gadget {
public:
	FOVAdjustmentGadget() :
		Gadget("FOVAdjustment", "gadgetron_mricore",
		"GenericReconFieldOfViewAdjustmentGadget")
	{
		add_property("debug_folder", "");
		add_property("perform_timing", "false");
		add_property("verbose", "false");
	}
};

class ScalingGadget : public Gadget {
public:
	ScalingGadget() :
		Gadget("Scaling", "gadgetron_mricore",
		"GenericReconImageArrayScalingGadget")
	{
		add_property("perform_timing", "false");
		add_property("verbose", "false");
		add_property("min_intensity_value", "64");
		add_property("max_intensity_value", "4095");
		add_property("scalingFactor", "10.0");
		add_property("use_constant_scalingFactor", "true");
		add_property("scalingFactor_dedicated", "100.0");
		add_property("auto_scaling_only_once", "true");
	}
};

class ImgArrSplitGadget : public Gadget {
public:
	ImgArrSplitGadget() :
		Gadget("ImageArraySplit", "gadgetron_mricore", "ImageArraySplitGadget")
	{}
};

class ExtGadget : public Gadget {
public:
	ExtGadget() :
		Gadget("Extract", "gadgetron_mricore", "ExtractGadget")
	{
		add_property("extract_mask", "1");
	}
};

class ComplexToFloatGadget : public Gadget {
public:
	ComplexToFloatGadget() :
		Gadget("ComplexToFloatAttrib", "gadgetron_mricore", "ComplexToFloatGadget")
	{}
};

class FloatToShortGadget : public Gadget {
public:
	FloatToShortGadget() :
		Gadget("FloatToShortAttrib", "gadgetron_mricore",
		"FloatToShortGadget")
	{
		add_property("min_intensity", "0");
		add_property("max_intensity", "32767");
		add_property("intensity_offset", "0");
	}
};

class ImgFinishGadget : public Gadget {
public:
	ImgFinishGadget() :
		Gadget("ImageFinish", "gadgetron_mricore", "ImageFinishGadget")
	{}
};

class AcqFinishGadget : public Gadget {
public:
	AcqFinishGadget() :
		Gadget("AcquisitionFinish", "gadgetron_mricore", "AcquisitionFinishGadget")
	{}
};

class SimpleReconstructionGadgetSet : public aGadget {
public:
	virtual void set_property(const char* prop, const char* value) 
	{
		if (boost::iequals(prop, "trigger_dimension") || 
			boost::iequals(prop, "sorting_dimension"))
			aat_.set_property(prop, value);
		else if (boost::iequals(prop, "n_dimension") || 
			boost::iequals(prop, "s_dimension") || 
			boost::iequals(prop, "split_slices"))
			bb_.set_property(prop, value);
		else
			THROW("unknown gadget parameter");
	}
	virtual std::string value_of(const char* prop)
	{
		if (boost::iequals(prop, "trigger_dimension") ||
			boost::iequals(prop, "sorting_dimension"))
			return aat_.value_of(prop);
		else if (boost::iequals(prop, "n_dimension") ||
			boost::iequals(prop, "s_dimension") ||
			boost::iequals(prop, "split_slices"))
			return bb_.value_of(prop);
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
};

#endif
