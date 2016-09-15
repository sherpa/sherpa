cp -r $RECIPE_DIR/../../* .

sed -i.orig "s|#install_dir=build|install_dir=$PREFIX|" setup.cfg


case $OSTYPE in
    darwin*)
#        export CFLAGS="-isysroot /Developer/SDKs/MacOSX10.5.sdk"

	# On Linux there is a libgfortran available in the environment, and statically linking it would require
	# rebuilding gcc from sources with PIC enabled. On OSX, no libgfortran is available, but statically
	# linking its dylib is not a problem.
        sed -i.orig "s|#extra-fortran-link-flags=|extra-fortran-link-flags=-static-libgfortran|" setup.cfg
	export LDFLAGS="-m64"
        cd extern
        ./configure --enable-standalone --enable-fftw --enable-region --enable-group --enable-wcs --enable-stk --prefix=$PREFIX --disable-maintainer-mode --enable-stuberrorlib --disable-shared --enable-shared=libgrp,stklib
        make
        make install
        touch built
        cd ..
        unset CFLAGS
        unset LDFLAGS
        ;;

    linux*)
        export CFLAGS="-L$PREFIX/lib"
        export LD_PRELOAD="$PREFIX/lib/libgfortran.so.3"
        ;;

esac

python setup.py clean --all

python setup.py install --prefix=$PREFIX

# This headers are known to collide with astropy's extensions and would prevent astropy from building
rm -f $PREFIX/include/wcs*.h

