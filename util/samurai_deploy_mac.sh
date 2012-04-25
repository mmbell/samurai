#!/bin/sh

rm -f lib/*
cp build/release/samurai bin/
cp -R /usr/local/Cellar/qt/4.7.4/lib/QtCore.framework/Versions/4/QtCore lib/QtCore
cp -R /usr/local/Cellar/qt/4.7.4/lib/QtGui.framework/Versions/4/QtGui lib/QtGui
cp -R /usr/local/Cellar/qt/4.7.4/lib/QtXml.framework/Versions/4/QtXml lib/QtXml
cp /usr/local/lib/libGeographic.dylib lib/libGeographic.dylib
cp /usr/local/lib/libhdf5.dylib lib/libhdf5.dylib
cp /usr/local/lib/libhdf5_hl.dylib lib/libhdf5_hl.dylib
cp /usr/local/lib/libnetcdf.dylib lib/libnetcdf.dylib
cp /usr/local/lib/libnetcdf_c++.dylib lib/libnetcdf_c++.dylib
cp /usr/lib/libcurl.4.dylib lib/libcurl.dylib
cp /usr/X11/lib/libpng15.15.dylib lib/libpng.dylib
cp /usr/local/lib/libsz.2.0.0.dylib lib/libsz.dylib

chmod 644 lib/*
install_name_tool -change /usr/local/lib/libhdf5.7.dylib @executable_path/../lib/libhdf5.dylib bin/samurai
install_name_tool -change /usr/local/lib/libnetcdf.7.dylib @executable_path/../lib/libnetcdf.dylib bin/samurai
install_name_tool -change /usr/local/lib/libnetcdf_c++.4.dylib @executable_path/../lib/libnetcdf_c++.dylib bin/samurai
install_name_tool -change /usr/local/lib/libGeographic.9.dylib @executable_path/../lib/libGeographic.dylib bin/samurai
install_name_tool -change /usr/lib/libcurl.4.dylib @executable_path/../lib/libcurl.dylib bin/samurai
install_name_tool -change /usr/local/Cellar/qt/4.7.4/lib/QtXml.framework/Versions/4/QtXml @executable_path/../lib/QtXml bin/samurai
install_name_tool -change /usr/local/Cellar/qt/4.7.4/lib/QtCore.framework/Versions/4/QtCore @executable_path/../lib/QtCore bin/samurai
install_name_tool -change /usr/local/Cellar/qt/4.7.4/lib/QtGui.framework/Versions/4/QtGui @executable_path/../lib/QtGui bin/samurai
install_name_tool -change /usr/X11/lib/libpng15.15.dylib @executable_path/../lib/libpng.dylib lib/QtGui
install_name_tool -change /usr/local/lib/libsz.2.0.0.dylib @executable_path/../lib/libsz.dylib lib/libhdf5.dylib
install_name_tool -change /usr/lib/libcurl.4.dylib @executable_path/../lib/libcurl.dylib lib/libnetcdf.dylib
install_name_tool -change /usr/local/lib/libsz.2.0.0.dylib @executable_path/../lib/libsz.dylib lib/libnetcdf.dylib
install_name_tool -change /usr/local/lib/libhdf5.7.dylib @executable_path/../lib/libhdf5.dylib lib/libnetcdf.dylib
install_name_tool -change /usr/local/lib/libhdf5_hl.7.dylib @executable_path/../lib/libhdf5_hl.dylib lib/libnetcdf.dylib
install_name_tool -change /usr/local/lib/libsz.2.0.0.dylib @executable_path/../lib/libsz.dylib lib/libhdf5_hl.dylib
install_name_tool -change /usr/local/Cellar/hdf5/1.8.8/lib/libhdf5.7.dylib @executable_path/../lib/libhdf5.dylib lib/libhdf5_hl.dylib
install_name_tool -change /usr/local/lib/libsz.2.0.0.dylib @executable_path/../lib/libsz.dylib lib/libnetcdf_c++.dylib
install_name_tool -change /usr/local/lib/libhdf5.7.dylib @executable_path/../lib/libhdf5.dylib lib/libnetcdf_c++.dylib
install_name_tool -change /usr/local/lib/libhdf5_hl.7.dylib @executable_path/../lib/libhdf5_hl.dylib lib/libnetcdf_c++.dylib
install_name_tool -change /usr/local/Cellar/netcdf/4.1.3/lib/libnetcdf.7.dylib @executable_path/../lib/libnetcdf.dylib lib/libnetcdf_c++.dylib
install_name_tool -change /usr/lib/libcurl.4.dylib @executable_path/../lib/libcurl.dylib lib/libnetcdf_c++.dylib

install_name_tool -id @executable_path/../lib/libhdf5.dylib ./lib/libhdf5.dylib
install_name_tool -id @executable_path/../lib/libhdf5_hl.dylib ./lib/libhdf5_hl.dylib
install_name_tool -id @executable_path/../lib/libnetcdf.dylib ./lib/libnetcdf.dylib
install_name_tool -id @executable_path/../lib/libnetcdf_c++.dylib ./lib/libnetcdf_c++.dylib
install_name_tool -id @executable_path/../lib/libGeographic.dylib ./lib/libGeographic.dylib
install_name_tool -id @executable_path/../lib/libcurl.dylib ./lib/libcurl.dylib
install_name_tool -id @executable_path/../lib/QtXml ./lib/QtXml
install_name_tool -id @executable_path/../lib/QtCore ./lib/QtCore
install_name_tool -id @executable_path/../lib/QtGui ./lib/QtGui
install_name_tool -id @executable_path/../lib/libpng.dylib ./lib/libpng.dylib
install_name_tool -id @executable_path/../lib/libsz.dylib ./lib/libsz.dylib

install_name_tool -change /usr/local/Cellar/qt/4.7.4/lib/QtCore.framework/Versions/4/QtCore @executable_path/../lib/QtCore ./lib/QtGui
install_name_tool -change /usr/local/Cellar/qt/4.7.4/lib/QtCore.framework/Versions/4/QtCore @executable_path/../lib/QtCore ./lib/QtXml

chmod 444 lib/*
rm samurai_v1.0.0-beta.tar.bz2
tar -cjvf samurai_v1.0.0-beta.tar.bz2 bin/ lib/ util/*.pl util/*.rb doc/ 

