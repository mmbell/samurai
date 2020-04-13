#ifndef STRINGS_H_
#define STRINGS_H_

#include <string>

static bool endsWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}


#endif // STRINGS_H_
