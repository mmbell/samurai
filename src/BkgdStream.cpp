#include <iostream>
#include "BkgdAdapter.h"

BkgdStream::BkgdStream(const char *fname)
{
  _stream = new std::ifstream(fname);
  if ( (_stream == NULL) || !_stream->good()) {
    std::cout << "Error opening " << fname << " for reading." << std::endl;
    exit(1);
  }
}

BkgdStream::~BkgdStream()
{
  delete _stream;
}

bool BkgdStream::next(int &time, real &lat, real &lon, real &alt, real &u,
		     real &v, real &w, real &t, real &qv, real &rhoa, real &qr)
{
  return *_stream >> time >> lat >> lon >> alt >> u >> v >> w >> t >> qv >> rhoa >> qr;
}
