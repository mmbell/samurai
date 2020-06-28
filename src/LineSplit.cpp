#include "LineSplit.h"
#include <sstream>
#include <string>

std::vector<std::string> LineSplit(const std::string& in, const char delim)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(in);
   while (std::getline(tokenStream, token, delim))
   {
      if (token != "") {
        tokens.push_back(token);
      }
   }
   return tokens;
}
