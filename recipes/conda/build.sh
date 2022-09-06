
sed -i.orig "s|#install_dir=build|install_dir=$PREFIX|" setup.cfg
git update-index --assume-unchanged setup.cfg

export PYTHON_LDFLAGS=" "

$PYTHON -m pip install --prefix=$PREFIX -vv .

# This headers are known to collide with astropy's extensions and would prevent astropy from building
rm -f $PREFIX/include/wcs*.h

