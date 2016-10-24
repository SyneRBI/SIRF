#ifndef GADGETS_LIBRARY
#define GADGETS_LIBRARY

#include <boost/algorithm/string.hpp>

#include "data_handle.h"

#define ADD_PROPERTY(P, V) \
	  xml_script += " <property>\n"; \
		xml_script += "  <name>" + std::string(P) + "</name>\n"; \
		xml_script += "  <value>" + V + "</value>\n"; \
		xml_script += " </property>\n";

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

class NoiseAdjustGadget : public aGadget {
public:
	NoiseAdjustGadget() {
		name_ = "NoiseAdjustGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
		xml_script += " <name>NoiseAdjust</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>NoiseAdjustGadget</classname>\n";
		xml_script += "</gadget>\n";
		return xml_script;
	}
};

class AsymmetricEchoGadget : public aGadget {
public:
	AsymmetricEchoGadget() {
		name_ = "AsymmetricEchoGadget";
	}
	virtual void set_property(const char* prop, const char* value) {}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
		xml_script += " <name>AsymmetricEcho</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>AsymmetricEchoAdjustROGadget</classname>\n";
		xml_script += "</gadget>\n";
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
	AcqAccTrigGadget() : 
		trigger_dimension_("repetition"), 
		sorting_dimension_("slice") 
	{
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
		ADD_PROPERTY("trigger_dimension", trigger_dimension_);
		ADD_PROPERTY("sorting_dimension", sorting_dimension_);
		xml_script += "</gadget>\n";
		return xml_script;
	}
private:
	std::string trigger_dimension_;
	std::string sorting_dimension_;
};

class BucketToBuffGadget : public aGadget {
public:
	BucketToBuffGadget() : 
		n_dimension_(""),
		s_dimension_(""), 
		split_slices_("true"),
		ignore_segment_("true"),
		verbose_("true")
	{
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
		else if (boost::iequals(prop, "ignore_segment"))
			ignore_segment_ = value;
		else if (boost::iequals(prop, "verbose"))
			verbose_ = value;
		else
			THROW("unknown gadget parameter");
	}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
    xml_script += " <name>Buff</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>BucketToBufferGadget</classname>\n";
		ADD_PROPERTY("verbose", verbose_);
		ADD_PROPERTY("N_dimension", n_dimension_);
		ADD_PROPERTY("S_dimension", s_dimension_);
		ADD_PROPERTY("split_slices", split_slices_);
		ADD_PROPERTY("ignore_segment", ignore_segment_);
		xml_script += "</gadget>\n";
		return xml_script;
	}
private:
	std::string n_dimension_;
	std::string s_dimension_;
	std::string split_slices_;
	std::string ignore_segment_;
	std::string verbose_;
};

class PrepRefGadget : public aGadget {
public:
	PrepRefGadget() : 
		debug_folder_(""),
		perform_timing_("true"),
		verbose_("true"),
		av_all_ref_N_("true"),
		av_all_ref_S_("true"),
		prep_ref_always_("true")
	{
		name_ = "PrepRefGadget";
	}
	virtual void set_property(const char* prop, const char* value)
	{
		if (boost::iequals(prop, "debug_folder"))
			debug_folder_ = value;
		else if (boost::iequals(prop, "perform_timing"))
			perform_timing_ = value;
		else if (boost::iequals(prop, "verbose"))
			verbose_ = value;
		else if (boost::iequals(prop, "average_all_ref_n"))
			av_all_ref_N_ = value;
		else if (boost::iequals(prop, "average_all_ref_s"))
			av_all_ref_S_ = value;
		else if (boost::iequals(prop, "prepare_ref_always"))
			prep_ref_always_ = value;
		else
			THROW("unknown gadget parameter");
	}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
		xml_script += " <name>PrepRef</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>GenericReconCartesianReferencePrepGadget</classname>\n";
		ADD_PROPERTY("debug_folder", debug_folder_);
		ADD_PROPERTY("perform_timing", perform_timing_);
		ADD_PROPERTY("verbose", verbose_);
		ADD_PROPERTY("average_all_ref_N", av_all_ref_N_);
		ADD_PROPERTY("average_all_ref_S", av_all_ref_S_);
		ADD_PROPERTY("prepare_ref_always", prep_ref_always_);
		xml_script += "</gadget>\n";
		return xml_script;
	}
private:
	std::string debug_folder_;
	std::string perform_timing_;
	std::string verbose_;
	std::string av_all_ref_N_;
	std::string av_all_ref_S_;
	std::string prep_ref_always_;
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

class CartesianGrappaGadget : public aGadget {
public:
	CartesianGrappaGadget() :
		image_series_("0"),
		coil_map_alg_("Inati"),
		dwnstr_coil_compr_("true"),
		dwnstr_coil_compr_th_("0.01"),
		dwnstr_coil_compr_nmod_("0"),
		debug_folder_(""),
		perform_timing_("true"),
		verbose_("true"),
		send_out_gfactor_("true")
	{
		name_ = "CartesianGrappaGadget";
	}
	virtual void set_property(const char* prop, const char* value)
	{
		if (boost::iequals(prop, "image_series"))
			image_series_ = value;
		else if (boost::iequals(prop, "coil_map_algorithm"))
			coil_map_alg_ = value;
		else if (boost::iequals(prop, "downstream_coil_compression"))
			dwnstr_coil_compr_ = value;
		else if (boost::iequals(prop, "downstream_coil_compression_thres"))
			dwnstr_coil_compr_th_ = value;
		else if (boost::iequals(prop, "downstream_coil_compression_modes"))
			dwnstr_coil_compr_nmod_ = value;
		else if (boost::iequals(prop, "debug_folder"))
			debug_folder_ = value;
		else if (boost::iequals(prop, "perform_timing"))
			perform_timing_ = value;
		else if (boost::iequals(prop, "verbose"))
			verbose_ = value;
		else if (boost::iequals(prop, "send_out_gfactor"))
			send_out_gfactor_ = value;
		else
			THROW("unknown gadget parameter");
	}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
		xml_script += " <name>CartesianGrappa</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>GenericReconCartesianGrappaGadget</classname>\n";
		ADD_PROPERTY("image_series", image_series_);
		ADD_PROPERTY("coil_map_algorithm", coil_map_alg_);
		ADD_PROPERTY("downstream_coil_compression", dwnstr_coil_compr_);
		ADD_PROPERTY("downstream_coil_compression_thres", dwnstr_coil_compr_th_);
		ADD_PROPERTY
			("downstream_coil_compression_num_modesKept", dwnstr_coil_compr_nmod_);
		ADD_PROPERTY("debug_folder", debug_folder_);
		ADD_PROPERTY("perform_timing", perform_timing_);
		ADD_PROPERTY("verbose", verbose_);
		ADD_PROPERTY("send_out_gfactor", send_out_gfactor_);
		xml_script += "</gadget>\n";
		return xml_script;
	}
private:
	std::string image_series_;
	std::string coil_map_alg_;
	std::string dwnstr_coil_compr_;
	std::string dwnstr_coil_compr_th_;
	std::string dwnstr_coil_compr_nmod_;
	std::string debug_folder_;
	std::string perform_timing_;
	std::string verbose_;
	std::string send_out_gfactor_;
};

class FOVAdjustmentGadget : public aGadget {
public:
	FOVAdjustmentGadget() :
		debug_folder_(""),
		perform_timing_("false"),
		verbose_("false")
	{
		name_ = "FOVAdjustmentGadget";
	}
	virtual void set_property(const char* prop, const char* value)
	{
		if (boost::iequals(prop, "debug_folder"))
			debug_folder_ = value;
		else if (boost::iequals(prop, "perform_timing"))
			perform_timing_ = value;
		else if (boost::iequals(prop, "verbose"))
			verbose_ = value;
		else
			THROW("unknown gadget parameter");
	}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
		xml_script += " <name>FOVAdjustment</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += 
			" <classname>GenericReconFieldOfViewAdjustmentGadget</classname>\n";
		ADD_PROPERTY("debug_folder", debug_folder_);
		ADD_PROPERTY("perform_timing", perform_timing_);
		ADD_PROPERTY("verbose", verbose_);
		xml_script += "</gadget>\n";
		return xml_script;
	}
private:
	std::string debug_folder_;
	std::string perform_timing_;
	std::string verbose_;
};

class ScalingGadget : public aGadget {
public:
	ScalingGadget() :
		perform_timing_("false"),
		verbose_("false"),
		min_intensity_value_("64"),
		max_intensity_value_("4095"),
		scalingFactor_("10.0"),
		use_constant_scalingFactor_("true"),
		scalingFactor_dedicated_("100.0"),
		auto_scaling_only_once_("true")
	{
		name_ = "ScalingGadget";
	}
	virtual void set_property(const char* prop, const char* value)
	{
		if (boost::iequals(prop, "perform_timing"))
			perform_timing_ = value;
		else if (boost::iequals(prop, "verbose"))
			verbose_ = value;
		else if (boost::iequals(prop, "min_intensity_value"))
			min_intensity_value_ = value;
		else if (boost::iequals(prop, "max_intensity_value"))
			max_intensity_value_ = value;
		else if (boost::iequals(prop, "scaling_factor"))
			scalingFactor_ = value;
		else if (boost::iequals(prop, "use_constant_scaling_factor"))
			use_constant_scalingFactor_ = value;
		else if (boost::iequals(prop, "scaling_factor_dedicated"))
			scalingFactor_dedicated_ = value;
		else if (boost::iequals(prop, "auto_scaling_only_once"))
			auto_scaling_only_once_ = value;
		else
			THROW("unknown gadget parameter");
	}
	virtual std::string xml() const {
		std::string xml_script("<gadget>\n");
		xml_script += " <name>Scaling</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script +=
			" <classname>GenericReconImageArrayScalingGadget</classname>\n";
		ADD_PROPERTY("perform_timing", perform_timing_);
		ADD_PROPERTY("verbose", verbose_);
		ADD_PROPERTY("min_intensity_value", min_intensity_value_);
		ADD_PROPERTY("max_intensity_value", max_intensity_value_);
		ADD_PROPERTY("scalingFactor", scalingFactor_);
		ADD_PROPERTY("use_constant_scalingFactor", use_constant_scalingFactor_);
		ADD_PROPERTY("scalingFactor_dedicated", scalingFactor_dedicated_);
		ADD_PROPERTY("auto_scaling_only_once", auto_scaling_only_once_);
		xml_script += "</gadget>\n";
		return xml_script;
	}
private:
	std::string perform_timing_;
	std::string verbose_;
	std::string min_intensity_value_;
	std::string max_intensity_value_;
	std::string scalingFactor_;
	std::string use_constant_scalingFactor_;
	std::string auto_scaling_only_once_;
	std::string scalingFactor_dedicated_;
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
	FloatToShortGadget() :
		min_intensity_("0"),
		max_intensity_("32767"),
		intensity_offset_("0")
	{
		name_ = "FloatToShortGadget";
	}
	virtual void set_property(const char* prop, const char* value) 
	{
		if (boost::iequals(prop, "min_intensity"))
			min_intensity_ = value;
		else if (boost::iequals(prop, "max_intensity"))
			max_intensity_ = value;
		else if (boost::iequals(prop, "intensity_offset"))
			intensity_offset_ = value;
		else
			THROW("unknown gadget parameter");
	}
	virtual std::string xml() const 
	{
		std::string xml_script("<gadget>\n");
		xml_script += " <name>FloatToShortAttrib</name>\n";
		xml_script += " <dll>gadgetron_mricore</dll>\n";
		xml_script += " <classname>FloatToUShortGadget</classname>\n";
		ADD_PROPERTY("min_intensity", min_intensity_);
		ADD_PROPERTY("max_intensity", max_intensity_);
		ADD_PROPERTY("intensity_offset", intensity_offset_);
		xml_script += "</gadget>\n";
		return xml_script;
	}
private:
	std::string min_intensity_;
	std::string max_intensity_;
	std::string intensity_offset_;
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
