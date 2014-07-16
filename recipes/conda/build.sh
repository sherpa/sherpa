cp -r $RECIPE_DIR/../../* .

sed -i.orig "s|#install_dir=build|install_dir=$PREFIX|" setup.cfg
sed -i.orig "s|#configure=None|configure=--disable-maintainer-mode --enable-shared --enable-stuberrorlib|" setup.cfg

export CFLAGS="-L$PREFIX/lib"

python setup.py install --prefix=$PREFIX
