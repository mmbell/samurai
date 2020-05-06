#include <string>
#include <unordered_map>

// Class to encapsulate error data read from Fractl

class ErrorData {

 public:
  
  ErrorData() { mishData = meshData = finalData = NULL; }
  ~ErrorData();
  
  double *init(std::string fname, std::unordered_map<std::string, std::string>* config, size_t numVar);

  double *getMeshVar(size_t var);	// After SA. Do we care anymore? TODO
  
  double *getMishData()  { return mishData; }
  double *getMeshData()  { return meshData; }  
  double *getFinalData() { return finalData; }
  
  void setMeshData(double *data)	{ meshData  = data; }
  void setFinalData(double *data)	{ finalData = data; }
  
  double meshValueAt(size_t va, size_t x, size_t y, size_t z);

  bool writeDebugNc(const std::string& netcdfFileName, bool mish, double *data);

  void getMishDims(size_t *dims) {	// TODO debug
    if (dims == NULL) return;
    dims[0] = mish_nz;
    dims[1] = mish_ny;
    dims[2] = mish_nx;
  }
  
 private:

  bool useDefaultValues;	// use values from the xml file, not from fractl
  size_t varDim;	// how many variables our caller thinks we have
  
  size_t mish_nx, mish_ny, mish_nz;	// size of the mish
  size_t mesh_nx, mesh_ny, mesh_nz;	// size of the mesh
  
  std::unordered_map <std::string, std::string> *configHash;
  double bgError[7];	// default fixed values from XML file.
  double *mishData;	// mish grid from the fractl ncd file
  double *meshData;	// DB + SA transformation of SBData
  double *finalData;	// std error version of finalAnalysis for Obs (after SI transform)
};
