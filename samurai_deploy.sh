#!/bin/sh

cp -R /Library/Frameworks/QtCore.framework/Versions/4/QtCore lib/QtCore
cp -R /Library/Frameworks/QtGui.framework/Versions/4/QtGui lib/QtGui
cp -R /Library/Frameworks/QtXml.framework/Versions/4/QtXml lib/QtXml
cp /Users/mmbell/Development/lib/libGeographic.5.dylib lib/libGeographic.5.dylib
cp /opt/local/lib/libcrypto.1.0.0.dylib lib/libcrypto.1.0.0.dylib
cp /opt/local/lib/libcurl.4.dylib lib/libcurl.4.dylib
cp /opt/local/lib/libhdf5.6.dylib lib/libhdf5.6.dylib
cp /opt/local/lib/libhdf5_hl.6.dylib lib/libhdf5_hl.6.dylib
cp /opt/local/lib/libiconv.2.dylib lib/libiconv.2.dylib
cp /opt/local/lib/libidn.11.dylib lib/libidn.11.dylib
cp /opt/local/lib/libintl.8.dylib lib/libintl.8.dylib
cp /opt/local/lib/libnetcdf.6.dylib lib/libnetcdf.6.dylib
cp /opt/local/lib/libnetcdf_c++.5.dylib lib/libnetcdf_c++.5.dylib
cp /opt/local/lib/libssl.1.0.0.dylib lib/libssl.1.0.0.dylib
cp /opt/local/lib/libsz.2.dylib lib/libsz.2.dylib
cp /opt/local/lib/libz.1.dylib lib/libz.1.dylib

install_name_tool -change /opt/local/lib/libcurl.4.dylib ./lib/libcurl.4.dylib samurai
install_name_tool -change /opt/local/lib/libhdf5.6.dylib ./lib/libhdf5.6.dylib samurai
install_name_tool -change /opt/local/lib/libnetcdf.6.dylib ./lib/libnetcdf.6.dylib samurai
install_name_tool -change /opt/local/lib/libnetcdf_c++.5.dylib ./lib/libnetcdf_c++.5.dylib samurai
install_name_tool -change /Users/mmbell/Development/lib/libGeographic.5.dylib ./lib/libGeographic.5.dylib samurai
install_name_tool -change QtXml.framework/Versions/4/QtXml ./lib/QtXml samurai
install_name_tool -change QtCore.framework/Versions/4/QtCore ./lib/QtCore samurai
install_name_tool -change QtGui.framework/Versions/4/QtGui ./lib/QtGui samurai

install_name_tool -id @executable_path/lib/libcurl.4.dylib ./lib/libcurl.4.dylib
install_name_tool -id @executable_path/lib/libhdf5.6.dylib ./lib/libhdf5.6.dylib
install_name_tool -id @executable_path/lib/libhdf5_hl.6.dylib ./lib/libhdf5_hl.6.dylib
install_name_tool -id @executable_path/lib/libnetcdf.6.dylib ./lib/libnetcdf.6.dylib
install_name_tool -id @executable_path/lib/libnetcdf_c++.5.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -id @executable_path/lib/libGeographic.5.dylib ./lib/libGeographic.5.dylib
install_name_tool -id @executable_path/lib/QtXml ./lib/QtXml
install_name_tool -id @executable_path/lib/QtCore ./lib/QtCore
install_name_tool -id @executable_path/lib/QtGui ./lib/QtGui
install_name_tool -id @executable_path/lib/libcrypto.1.0.0.dylib ./lib/libcrypto.1.0.0.dylib
install_name_tool -id @executable_path/lib/libidn.11.dylib ./lib/libidn.11.dylib
install_name_tool -id @executable_path/lib/libintl.8.dylib ./lib/libintl.8.dylib
install_name_tool -id @executable_path/lib/libiconv.2.dylib ./lib/libiconv.2.dylib
install_name_tool -id @executable_path/lib/libssl.1.0.0.dylib ./lib/libssl.1.0.0.dylib
install_name_tool -id @executable_path/lib/libsz.2.dylib ./lib/libsz.2.dylib
install_name_tool -id @executable_path/lib/libsz.1.dylib ./lib/libz.1.dylib

install_name_tool -change QtCore.framework/Versions/4/QtCore @executable_path/lib/QtCore ./lib/QtGui
install_name_tool -change QtCore.framework/Versions/4/QtCore @executable_path/lib/QtCore ./lib/QtXml

install_name_tool -change /opt/local/lib/libintl.8.dylib @executable_path/lib/libintl.8.dylib ./lib/libidn.11.dylib
install_name_tool -change /opt/local/lib/libiconv.2.dylib @executable_path/lib/libiconv.2.dylib ./lib/libidn.11.dylib
install_name_tool -change /opt/local/lib/libiconv.2.dylib @executable_path/lib/libiconv.2.dylib ./lib/libintl.8.dylib

install_name_tool -change /opt/local/lib/libidn.11.dylib @executable_path/lib/libidn.11.dylib ./lib/libcurl.4.dylib
install_name_tool -change /opt/local/lib/libcrypto.1.0.0.dylib @executable_path/lib/libcrypto.1.0.0.dylib  ./lib/libcurl.4.dylib
install_name_tool -change /opt/local/lib/libssl.1.0.0.dylib @executable_path/lib/libssl.1.0.0.dylib ./lib/libcurl.4.dylib
install_name_tool -change /opt/local/lib/libz.1.dylib @executable_path/lib/libz.1.dylib ./lib/libcurl.4.dylib

install_name_tool -change /opt/local/lib/libz.1.dylib @executable_path/lib/libz.1.dylib ./lib/libcrypto.1.0.0.dylib
install_name_tool -change /opt/local/lib/libcurl.4.dylib @executable_path/lib/libcurl.4.dylib ./lib/libnetcdf.6.dylib
install_name_tool -change /opt/local/lib/libidn.11.dylib @executable_path/lib/libidn.11.dylib ./lib/libnetcdf.6.dylib
install_name_tool -change /opt/local/lib/libintl.8.dylib @executable_path/lib/libintl.8.dylib ./lib/libnetcdf.6.dylib
install_name_tool -change /opt/local/lib/libiconv.2.dylib @executable_path/lib/libiconv.2.dylib ./lib/libnetcdf.6.dylib
install_name_tool -change /opt/local/lib/libssl.1.0.0.dylib @executable_path/lib/libssl.1.0.0.dylib ./lib/libnetcdf.6.dylib
install_name_tool -change /opt/local/lib/libcrypto.1.0.0.dylib @executable_path/lib/libcrypto.1.0.0.dylib ./lib/libnetcdf.6.dylib
install_name_tool -change /opt/local/lib/libhdf5_hl.6.dylib @executable_path/lib/libhdf5_hl.6.dylib ./lib/libnetcdf.6.dylib
install_name_tool -change /opt/local/lib/libhdf5.6.dylib @executable_path/lib/libhdf5.6.dylib ./lib/libnetcdf.6.dylib
install_name_tool -change /opt/local/lib/libsz.2.dylib @executable_path/lib/libsz.2.dylib ./lib/libnetcdf.6.dylib
install_name_tool -change /opt/local/lib/libz.1.dylib @executable_path/lib/libz.1.dylib ./lib/libnetcdf.6.dylib

install_name_tool -change /opt/local/lib/libnetcdf.6.dylib @executable_path/lib/libnetcdf.6.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libcurl.4.dylib @executable_path/lib/libcurl.4.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libidn.11.dylib @executable_path/lib/libidn.11.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libintl.8.dylib @executable_path/lib/libintl.8.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libiconv.2.dylib @executable_path/lib/libiconv.2.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libssl.1.0.0.dylib @executable_path/lib/libssl.1.0.0.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libcrypto.1.0.0.dylib @executable_path/lib/libcrypto.1.0.0.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libhdf5_hl.6.dylib @executable_path/lib/libhdf5_hl.6.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libhdf5.6.dylib @executable_path/lib/libhdf5.6.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libsz.2.dylib @executable_path/lib/libsz.2.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libz.1.dylib @executable_path/lib/libz.1.dylib ./lib/libnetcdf_c++.5.dylib

install_name_tool -change /opt/local/lib/libsz.2.dylib @executable_path/lib/libsz.2.dylib ./lib/libhdf5.6.dylib
install_name_tool -change /opt/local/lib/libz.1.dylib @executable_path/lib/libz.1.dylib ./lib/libhdf5.6.dylib

install_name_tool -change /opt/local/lib/libsz.2.dylib @executable_path/lib/libsz.2.dylib ./lib/libhdf5_hl.6.dylib
install_name_tool -change /opt/local/lib/libz.1.dylib @executable_path/lib/libz.1.dylib ./lib/libhdf5_hl.6.dylib

install_name_tool -change /opt/local/lib/libidn.11.dylib @executable_path/lib/libidn.11.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libiconv.2.dylib @executable_path/lib/libiconv.2.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libssl.1.0.0.dylib @executable_path/lib/libssl.1.0.0.dylib ./lib/libnetcdf_c++.5.dylib
install_name_tool -change /opt/local/lib/libcrypto.1.0.0.dylib @executable_path/lib/libcrypto.1.0.0.dylib ./lib/libnetcdf_c++.5.dylib

install_name_tool -change /opt/local/lib/libz.1.dylib @executable_path/lib/libz.1.dylib ./lib/libssl.1.0.0.dylib
install_name_tool -change /opt/local/lib/libcrypto.1.0.0.dylib @executable_path/lib/libcrypto.1.0.0.dylib ./lib/libssl.1.0.0.dylib

rm samurai_v0.2.1.tar.bz2
tar -cjvf samurai_v0.2.1.tar.bz2 ecgrads2bg.pl met_formulas.pl samurai samurai_Background.in samurai_convective.xml lib/ vardata/

