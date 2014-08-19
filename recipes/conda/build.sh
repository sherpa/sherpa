cp -r $RECIPE_DIR/../../* .

sed -i.orig "s|#install_dir=build|install_dir=$PREFIX|" setup.cfg

case $OSTYPE in
    darwin*)
        export CFLAGS="-isysroot /Developer/SDKs/MacOSX10.5.sdk"
        export LDFLAGS="-m64"
        cd extern
        ./configure --enable-fftw --enable-region --enable-group --enable-wcs --prefix=$PREFIX --disable-maintainer-mode --enable-stuberrorlib --disable-shared --enable-shared=libgrp
        make
        make install
        touch built
        cd ..
        export LDFLAGS="-undefined dynamic_lookup"
        ;;

    linux*)
        export CFLAGS="-L$PREFIX/lib"
        ;;

esac

unset CFLAGS
unset LDFLAGS

python setup.py install --prefix=$PREFIX
