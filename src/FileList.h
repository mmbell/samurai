#ifndef FILELIST_H_
#define FILELIST_H_

#include <vector>
#include <string>

std::vector<std::string> FileList(std::string &path_in);

// Get a file's extension, since we can't use std::filesystem w/ Intel c++
std::string Extension(std::string);

// Check if a directory exists:
bool DirectoryExists(std::string);

#endif // FILELIST_H_
