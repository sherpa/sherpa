install-exec-local: group.so
	cp -f group.so $(libdir)/

uninstall-local:
	rm -f $(libdir)/group.so

group.so : 
	python setup.py build_ext --inplace

clean-local:
	rm -rf group.so build

