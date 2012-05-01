rm -f lib/*
cp build/release/samurai bin/
cp /work/mmbell/lib/libQtCore.so.4 lib/
cp /work/mmbell/lib/libQtGui.so.4 lib/
cp /work/mmbell/lib/libQtXml.so.4 lib/
cp /work/mmbell/lib/libcurl.so.4 lib/
cp /work/mmbell/lib/libhdf5.so.6 lib/
cp /work/mmbell/lib/libnetcdf.so.7 lib/
cp /work/mmbell/lib/libnetcdf_c++.so.4 lib/
cp /work/mmbell/lib/libGeographic.so.5 lib/
cp /opt/intel/composerxe-2011.4.191/compiler/lib/intel64/libiomp5.so lib/
cp /work/mmbell/lib/libhdf5_hl.so.6 lib/
chmod 444 lib/*
rm samurai_v1.0.0-RC1_Linux_intel.tar.bz2
tar -cjvf samurai_v1.0.0-RC1_Linux_intel.tar.bz2 bin/ lib/ util/ecgrads2bg.pl util/samurai_out2bg.pl util/samurai_lineartrack.pl doc/ 

