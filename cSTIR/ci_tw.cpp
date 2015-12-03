#include <fstream>
#include <string.h>

#include "TextWriter.h"

extern "C" {
	void* newTextPrinter(const char* stream) {
		TextPrinter* p = new TextPrinter(stream);
		return (void*)p;
	}
	void* newTextWriter(const char* stream) {
		TextWriter* w = new TextWriter;
		w->out = new std::ofstream;
		((std::ofstream*)w->out)->open(stream, std::ios::out);
		return (void*)w;
	}
	void openChannel(int channel, void* ptr_w) {
		TextWriterHandle h;
		switch (channel) {
		case INFORMATION_CHANNEL:
			h.set_information_channel((aTextWriter*)ptr_w);
			break;
		case WARNING_CHANNEL:
			h.set_warning_channel((aTextWriter*)ptr_w);
			break;
		case ERROR_CHANNEL:
			h.set_error_channel((aTextWriter*)ptr_w);
			break;
		default:
			h.set_information_channel((aTextWriter*)ptr_w);
			h.set_warning_channel((aTextWriter*)ptr_w);
			h.set_error_channel((aTextWriter*)ptr_w);
		}
	}
	void closeChannel(int channel, void* ptr_w) {
		TextWriterHandle h;
		switch (channel) {
		case INFORMATION_CHANNEL:
			if (h.information_channel_ptr() == ptr_w)
				h.set_information_channel(0);
			break;
		case WARNING_CHANNEL:
			if (h.warning_channel_ptr() == ptr_w)
				h.set_warning_channel(0);
			break;
		case ERROR_CHANNEL:
			if (h.error_channel_ptr() == ptr_w)
				h.set_error_channel(0);
			break;
		default:
			h.set_information_channel(0);
			h.set_warning_channel(0);
			h.set_error_channel(0);
		}
	}
	void setWriter(void* ptr_w, int channel) {
		TextWriterHandle h;
		switch (channel) {
		case INFORMATION_CHANNEL:
			h.set_information_channel((aTextWriter*)ptr_w);
			break;
		case WARNING_CHANNEL:
			h.set_warning_channel((aTextWriter*)ptr_w);
			break;
		case ERROR_CHANNEL:
			h.set_error_channel((aTextWriter*)ptr_w);
			break;
		}
	}
	void resetWriter() {
		TextWriterHandle h;
		h.set_information_channel(0);
		h.set_warning_channel(0);
		h.set_error_channel(0);
	}
	void deleteTextPrinter(void* ptr) {
		delete (TextPrinter*)ptr;
	}
	void deleteTextWriter(void* ptr_w) {
		if (!ptr_w)
			return;
		TextWriter* w = (TextWriter*)ptr_w;
		if (!w->out)
			return;
		((std::ofstream*)w->out)->close();
		((std::ofstream*)w->out)->clear();
		delete w->out;
		delete w;
	}
}
