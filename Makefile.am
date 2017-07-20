## Helloworld C++: an example project using Automake

## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}
AUTOMAKE_OPTIONS = foreign
SUBDIRS = src data scripts

include Makefile.vars

install-exec-hook:
	cd data && rm -rf $(CONFIG) && tar xzvf $(CONFIG).tar.gz; cd $(CONFIG) && perl Makefile.PL INSTALL_BASE=$(prefix) && make && make test && make install && cd .. && rm -rf $(CONFIG)
