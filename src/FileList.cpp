#include "FileList.h"

#include <iostream>
#include <algorithm>

#include <dirent.h> // Linux tool, not standard C++

std::vector<std::string> FileList(std::string &path_in)
{
	// Set up a vector of paths:
	std::vector<std::string> list;

	// Traverse the directory:
  auto directory = opendir(path_in.c_str());

  struct dirent *dp;
  while ((dp = readdir(directory)) != NULL) {
     //std::cout << "Found file : " << dp->d_name << std::endl;
     if (dp->d_type == DT_REG) {
      list.push_back(dp->d_name);
    }
  }


	// Sort it:
	std::sort(list.begin(), list.end());

	// Debug
	std::cout << "Debug: num files = " << list.size() << std::endl;

  closedir(directory);

	// return the vector
	return list;
}


std::string Extension(std::string in) {
  auto index = in.rfind('.');
  //std::cout << "Doing extension of: " << in << " ( index = " << index << " ) " << std::endl;
  if (index != std::string::npos) {
    return in.substr(index);
  } 
  return "";
}


bool DirectoryExists(std::string in) {
   auto directory = opendir(in.c_str());

   if (directory != NULL) {
			closedir(directory);
      return true;
   }

	closedir(directory);

   return false;
}

