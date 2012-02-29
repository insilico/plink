ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign
AM_CXXFLAGS = @OS_CXX@
lib_LTLIBRARIES=libplink.la
libplink_la_SOURCES=\
options.cpp input.cpp binput.cpp tinput.cpp genome.cpp	\
helper.cpp stats.cpp filters.cpp locus.cpp multi.cpp crandom.cpp	\
cluster.cpp mds.cpp output.cpp informative.cpp assoc.cpp epi.cpp	\
prephap.cpp phase.cpp trio.cpp tdt.cpp sharing.cpp genepi.cpp sets.cpp	\
perm.cpp mh.cpp genedrop.cpp gxe.cpp merge.cpp hotel.cpp multiple.cpp	\
haploCC.cpp haploTDT.cpp poo.cpp webcheck.cpp qfam.cpp linear.cpp	\
bmerge.cpp parse.cpp mishap.cpp legacy.cpp homozyg.cpp segment.cpp	\
model.cpp logistic.cpp glm.cpp dcdflib.cpp elf.cpp dfam.cpp fisher.cpp	\
linput.cpp sockets.cpp lookup.cpp proxy.cpp pdriver.cpp haploQTL.cpp	\
haplohelper.cpp haplowindow.cpp genogroup.cpp nonfounderphasing.cpp	\
clumpld.cpp genoerr.cpp em.cpp impute.cpp metaem.cpp profile.cpp	\
nlist.cpp whap.cpp simul.cpp gvar.cpp cnv.cpp step.cpp greport.cpp	\
flip.cpp qualscores.cpp cnvqt.cpp cfamily.cpp setscreen.cpp idhelp.cpp	\
tag.cpp hapglm.cpp lookup2.cpp blox.cpp zed.cpp dosage.cpp annot.cpp	\
metaanal.cpp zfstream.cpp plinklibhandler.cpp

libplink_la_HEADERS=\
plink.h options.h helper.h stats.h crandom.h sets.h phase.h	\
perm.h model.h linear.h logistic.h dcdflib.h ipmpar.h cdflib.h		\
fisher.h sockets.h haplowindow.h genogroup.h clumpld.h nlist.h whap.h	\
gvar.h cnv.h cfamily.h idhelp.h zed.h zfstream.h plinklibhandler.h

libplink_la_LDFLAGS = -release 1.0.7

if WITH_R_PLUGINS
libplink_la_SOURCES += r.cpp Rconnection.cpp
libplink_la_HEADERS += sisocks.h Rsrv.h Rconnection.h config.h
AM_CXXFLAGS += -DWITH_R_PLUGINS
LIBS += -ldl
if SYS_UNIX
LIBS += -lcrypt
endif
endif

if WITH_WEBCHECK
if SYS_WIN
LIBS += -lwsock32
endif
else
AM_CXXFLAGS += -DSKIP
endif

if WITH_LAPACK
libplink_la_SOURCES += lapackf.cpp
libplink_la_HEADERS += lapackf.h
AM_CXXFLAGS += -DWITH_LAPACK
endif

if WITH_OPENMP
AM_CXXFLAGS += -fopenmp
endif

if FORCE_32BIT
AM_CXXFLAGS += -m32
endif

if FORCE_DYNAMIC

else
AM_CXXFLAGS += -static
endif

# always build statically on Windows
if SYS_WIN
AM_CXXFLAGS += -static -static-libgcc -static-libstdc++
endif

# OS X flags
if SYS_MAC
AM_CXXFLAGS += -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen
endif

# dir to install headers for plink library
libplink_ladir=$(includedir)/plink

# build plink binary, linking to library built above
bin_PROGRAMS = plink
plink_SOURCES =  plink.cpp
plink_LDADD= -lplink
