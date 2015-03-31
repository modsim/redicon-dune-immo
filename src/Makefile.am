SUBDIRS = tests

noinst_PROGRAMS = immoplate immoplate-bridge
EXTRA_DIST = immo-run immo-analyze.sh

immoplate_SOURCES = immoplate.cc boundary-robin.hh boundary-dirichelt.hh\
		lib/assembleMatrix.hh lib/updateMatrix.hh lib/assembleVector.hh \
		lib/saveField.hh lib/saveField2.hh lib/index.hh\
		lib/laplacian.hh lib/gamma.hh lib/sourceFunctions.hh\
		lib/integrateField.hh lib/integrate_n.hh lib/integrateEntity.hh \
		lib/evalField.hh lib/calcSurfArea.hh\
		lib/setDirichlet.hh lib/setRobin.hh \
		lib/integrateEntityDynamic.hh \
		lib/pointSource.hh lib/diracDelta.hh lib/insideEntity.hh \
		lib/entityFunction.hh lib/entityFunction2.hh lib/shapeFunctions.hh 

immoplate_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS) \
	-I $(top_srcdir)/src/lib/

# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
immoplate_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) \
	$(DUNEMPILIBS)	\
	$(LDADD)
immoplate_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)

immoplate_bridge_SOURCES = immoplate-bridge.cc bridge-forces.hh \
		 bridge.hh bridge-entityFunction.hh bridge-robinBC.hh \
		 boundary-robin.hh boundary-dirichelt.hh\
		lib/assembleMatrix.hh lib/updateMatrix.hh lib/assembleVector.hh \
		lib/saveField.hh lib/saveField2.hh lib/index.hh\
		lib/laplacian.hh lib/gamma.hh lib/sourceFunctions.hh\
		lib/integrateField.hh lib/integrate_n.hh lib/integrateEntity.hh \
		lib/evalField.hh lib/calcSurfArea.hh\
		lib/setDirichlet.hh lib/setRobin.hh \
		lib/integrateEntityDynamic.hh \
		lib/pointSource.hh lib/diracDelta.hh lib/insideEntity.hh \
		lib/entityFunction.hh lib/entityFunction2.hh lib/shapeFunctions.hh 

immoplate_bridge_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS) \
	-I $(top_srcdir)/src/lib/
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
immoplate_bridge_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) \
	$(DUNEMPILIBS)	\
	$(LDADD)
immoplate_bridge_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)


# don't follow the full GNU-standard
# we need automake 1.9
AUTOMAKE_OPTIONS = foreign 1.9

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-geometry=$(DUNE_GEOMETRY_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT)  CXX="$(CXX)" CC="$(CC)"

EXTRA_DIST = CMakeLists.txt

include $(top_srcdir)/am/global-rules

