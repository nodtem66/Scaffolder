#include "utils.h"

int util::make_dir(std::string& str) {
    if (str.empty())
        return 0;
#ifdef _WIN32
    return _mkdir(str.c_str());
#else
    return mkdir(str.c_str(), 0733);
#endif
}

    /* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
     * windows and linux. */

util::uint64 util::GetTimeMs64()
{
#ifdef _WIN32
    /* Windows */
    FILETIME ft;
    LARGE_INTEGER li{};

    /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
        * to a LARGE_INTEGER structure. */
    GetSystemTimeAsFileTime(&ft);
    li.LowPart = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;

    util::uint64 ret = li.QuadPart;
    ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
    ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

    return ret;
#else
    /* Linux */
    struct timeval tv;

    gettimeofday(&tv, NULL);

    uint64 ret = tv.tv_usec;
    /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
    ret /= 1000;

    /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
    ret += (tv.tv_sec * 1000);

    return ret;
#endif
 }

std::string util::PathGetBasename(std::string const& path)
{
    size_t firstindex = path.find_last_of("/\\");
    firstindex = firstindex == std::string::npos ? 0 : firstindex + 1;
    size_t lastindex = path.find_last_of(".");
    return path.substr(firstindex, lastindex - firstindex);
}

std::string util::PathGetExtension(std::string const& path)
{
    std::string ext;

    // Find the last dot, if any.
    size_t dotIdx = path.find_last_of(".");
    if (dotIdx != std::string::npos)
    {
        // Find the last directory separator, if any.
        size_t dirSepIdx = path.find_last_of("/\\");

        // If the dot is at the beginning of the file name, do not treat it as a file extension.
        // e.g., a hidden file:  ".alpha".
        // This test also incidentally avoids a dot that is really a current directory indicator.
        // e.g.:  "alpha/./bravo"
        if (dotIdx > dirSepIdx + 1)
        {
            ext = path.substr(dotIdx);
        }
    }

    return ext;
}

void util::to_lower(std::string& s) {
    std::locale loc;
    for (std::string::size_type i = 0; i < s.length(); ++i)
        s[i] = std::tolower(s[i], loc);
}