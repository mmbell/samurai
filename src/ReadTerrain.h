/*
 *  ReadTerrain.h
 *  samurai
 *
 *  Created by Ting-Yu Cha on 1/16/2022.
 *  Copyright 2022 Ting-Yu Cha. All rights reserved.
 *
 */

// #ifndef READTERRAIN_H
// #define READTERRAIN_H

#include "VarDriver.h"
#include "VarDriver3D.h"
#include "MetObs.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

class ReadTerrain : public VarDriver3D
{
public:
    bool readTerrainTXT(std::string &filename, std::vector<MetObs>* terrainData)
  {
    std::ifstream metFile(filename);
    if (!metFile.is_open()) {
      return false;
  	}
    MetObs ob;
    std::string line2;

    while (std::getline(metFile, line2)) {
  		// ob.setTime(datetime_);
      auto parts = LineSplit(line2, ' ');
      ob.setLat(std::stof(parts[0]));
      ob.setLon(std::stof(parts[1]));
      ob.setAltitude(std::stof(parts[2]));
      ob.setTerrainDX(std::stof(parts[3]));
      ob.setTerrainDY(std::stof(parts[4]));
      ob.setTerrainX(std::stof(parts[5]));
      ob.setTerrainY(std::stof(parts[6]));
      ob.setObType(MetObs::terrain);
      terrainData->push_back(ob);
    }
  	// std::cout << "Successfully read the terrain file" << std::endl;
    metFile.close();
    return true;
    // return read_terrain(metFile, metData);
  }
};
