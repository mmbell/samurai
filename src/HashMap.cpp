
#include "HashMap.h"

HashMap::HashMap(void)
{

}

std::string HashMap::operator [](const std::string& key)
{
  // This function duplicates Qt's HashMap operations by returning a "0" value when a key is not found:
  auto value = map_.find(key);

  if (value == map_.end()) {
    return "0";
  }

  return value->second;
}

bool HashMap::insert(const std::string& key, const std::string& value)
{
  map_.insert({key, value});
  return true;
}

bool HashMap::exists(const std::string& key)
{
  if (map_.find(key) != map_.end()) {
    return true;
  } 
  return false;
}

