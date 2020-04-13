#include <iostream>
#include "Xml.h"


std::vector<const XMLElement* > XMLGetElements(const XMLNode* config)
{
	std::vector<const XMLElement* > list;

	for (auto elem = config->FirstChildElement(); elem != NULL; elem = elem->NextSiblingElement()) {
		list.push_back(elem);
  }

	return list;
}

std::vector<const XMLAttribute* > XMLGetAttributes(const XMLElement* element)
{
	std::vector<const XMLAttribute* > list;

  for (auto attribute = element->FirstAttribute(); attribute != NULL; attribute = attribute->Next()) {
		list.push_back(attribute);
  }

	return list;
}

XMLDocument readXmlConfig(const char *path)
{
  XMLDocument foo;
}

