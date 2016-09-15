/usr/bin/perl /usr/share/perl/5.10/ExtUtils/xsubpp  -typemap /usr/share/perl/5.10/ExtUtils/typemap   hypergeom_bf98.xs > hypergeom_bf98.xsc && mv hypergeom_bf98.xsc hypergeom_bf98.c
cc -c  -I/home/pignatelli/workspace/CROC-standAlone -D_REENTRANT -D_GNU_SOURCE -DDEBIAN -fno-strict-aliasing -pipe -I/usr/local/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O2 -g   -DVERSION=\"0.00\" -DXS_VERSION=\"0.00\" -fPIC "-I/usr/lib/perl/5.10/CORE"   hypergeom_bf98.c
hypergeom_bf98.xs:7:1: warning: "MIN" redefined
In file included from /usr/lib/perl/5.10/CORE/perl.h:653,
                 from hypergeom_bf98.xs:2:
/usr/include/sys/param.h:74:1: warning: this is the location of the previous definition
hypergeom_bf98.xs:8:1: warning: "MAX" redefined
/usr/include/sys/param.h:75:1: warning: this is the location of the previous definition
Running Mkbootstrap for hypergeom_bf98 ()
chmod 644 hypergeom_bf98.bs
rm -f blib/arch/auto/hypergeom_bf98/hypergeom_bf98.so
cc  -shared -O2 -g -L/usr/local/lib hypergeom_bf98.o  -o blib/arch/auto/hypergeom_bf98/hypergeom_bf98.so 	\
	     	\
	  
chmod 755 blib/arch/auto/hypergeom_bf98/hypergeom_bf98.so
cp hypergeom_bf98.bs blib/arch/auto/hypergeom_bf98/hypergeom_bf98.bs
chmod 644 blib/arch/auto/hypergeom_bf98/hypergeom_bf98.bs
