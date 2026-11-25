#include "CleanFiles.h"
#include <stdexcept>
bool clear_directory(const std::string&) { throw std::runtime_error("Clearing output is disabled; create a fresh run directory."); }
