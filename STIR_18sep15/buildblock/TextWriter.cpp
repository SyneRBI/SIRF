#include "TextWriter.h"

aTextWriter* TextWriterHandle::writer = 0;
char Stir::endl = '\n';
AltCout Stir::cout;
AltCout Stir::cerr;

void writeText(const char* text) {
	TextWriterHandle h;
	h.write(text);
}

