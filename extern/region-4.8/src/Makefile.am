BUILT_SOURCES = regparser.c reglexer.c

# Overrides the prefix name in generated bison code. Replace if upgrading to bison 2.6+
AM_YFLAGS = -d -p"regYY"

# Need this for ylwrap
AM_LFLAGS = -o lex.yy.c

lib_LTLIBRARIES = libregion.la

include_HEADERS = cxcregion.h region_priv.h regpar.h

libregion_la_SOURCES = \
    reglexer.l regparser.y \
 	region_area.c \
	region_convert.c \
	region_create.c \
	region_extent.c \
	region_eval.c \
 	region_logic.c \
	region_repr.c \
	region_shape.c \
 	regShape_utils.c \
 	regAnnulus.c regField.c regPolygon.c regSector.c regMask.c \
 	regCircle.c regEllipse.c regBox.c regPoint.c regRectangle.c regPie.c

