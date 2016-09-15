cp -r $RECIPE_DIR/../../* .

sed -i.orig "s|#install_dir=build|install_dir=$PREFIX|" setup.cfg

python setup.py clean --all

case $OSTYPE in
    darwin*)
    export LDFLAGS="-m64"
    ;;

    linux*)
    export CFLAGS="-L$PREFIX/lib"
    ;;

esac

python setup.py install --prefix=$PREFIX

# This headers are known to collide with astropy's extensions and would prevent astropy from building
rm -f $PREFIX/include/wcs*.h

