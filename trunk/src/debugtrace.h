#ifndef DEBUGTRACE_H
#define DEBUGTRACE_H
#include <cstdio>
#include <exception>
#include <cstdarg>
namespace meshtalent {
	class DebugTrace {
	public:
		DebugTrace(char* tracefilename, bool isAppend = false) : ptracefile(NULL) {
			if (isAppend) {
				ptracefile = fopen(tracefilename, "a");
			} else {
				ptracefile = fopen(tracefilename, "w");
			}
			if (!ptracefile) {
				throw std::exception();
			}
		}
		~DebugTrace() {
			fclose(ptracefile);
		}
		int Trace(const char* format, ...) {
			va_list ap;
			int n;
			va_start(ap, format);
			n = vfprintf(ptracefile, format, ap);
			va_end(ap);
			return n;
		}
	private:
		FILE* ptracefile;
	};
} // end of namespace meshtalent

#endif // DEBUGTRACE_H
