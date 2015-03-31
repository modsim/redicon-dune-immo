check_PROGRAMS = test-bridge  test-bridge-force

INCLUDES = -I$(top_srcdir)/src/lib/ -I$(top_srcdir)/src

AM_CPPFLAGS= -Wall -DMAINTAINER=\"@MAINTAINER@\"
AM_LDFLAGS =  -Wall -pedantic -ansi

test_bridge_SOURCES = testBridge.cc 
test_bridge_CPPFLAGS = $(AM_CPPFLAGS) \
        $(DUNEMPICPPFLAGS) \
        $(UG_CPPFLAGS) \
        $(AMIRAMESH_CPPFLAGS) \
        $(ALBERTA_CPPFLAGS) \
        $(ALUGRID_CPPFLAGS) \
        -I $(top_srcdir)/src/lib/

test_bridge_force_SOURCES = testBridgeForce.cc 
test_bridge_force_CPPFLAGS = $(AM_CPPFLAGS) \
        $(DUNEMPICPPFLAGS) \
        $(UG_CPPFLAGS) \
        $(AMIRAMESH_CPPFLAGS) \
        $(ALBERTA_CPPFLAGS) \
        $(ALUGRID_CPPFLAGS) \
        -I $(top_srcdir)/src/lib/

