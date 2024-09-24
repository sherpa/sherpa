export PYTHON_LDFLAGS=" "

$PYTHON -m pip install --prefix=$PREFIX -vv .

# This headers are known to collide with astropy's extensions and would prevent astropy from building
rm -f $PREFIX/include/wcs*.h

