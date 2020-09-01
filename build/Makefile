
ifndef $(MACH)
	#MACH=MAC
	MACH=FSL
endif

include make.$(MACH)

####################   DIRECTORIES   ##########################################

SRCDIR       = ../source/
CVODE_SRCDIR = $(SRCDIR)cvode/
PROG         = ../run/odt.x

YAML_LIB = ../source/yaml/lib/
YAML_INC = ../source/yaml/include/

####################   COMPILATION   ##########################################

includes += -I$(SRCDIR) -I$(SRCDIR)flmlt -I$(SRCDIR)hips -I$(SRCDIR)linevariables/
includes += -I$(SRCDIR)linevariables/soot/
includes += -I$(SRCDIR)odtcases/ -I$(SRCDIR)nr/

####################   Cantera flags   ########################################

ext_libs += -lcantera -lpthread
includes += -I$(CANT_INC)
libs     += -L$(CANT_LIB)

####################   YAML flags   ########################################

ext_libs += -lyaml-cpp -lboost_system -lblas -llapack
includes += -I$(YAML_INC) -I/usr/local/include
libs     += -L$(YAML_LIB)

####################   other flags   ##########################################
ifneq ($(CHEMISTRY), CANTERARR)
    D1        = -DPROBLEMSPECIFICRR
endif
ifeq ($(PARALLEL), YES)
    D2        = -DDOMPI
    includes += -I$(MPI_PATH)include/
endif
ifeq ($(SILENT), YES)
    D3        = -DSILENT
endif

DIRECTIVES = $(D1) $(D2) $(D3) -D_USESTDVECTOR_


COMPCC   = $(CC_COMPILER_NAME) $(includes) $(options) $(cppoptions) $(DIRECTIVES) -c
COMPC    = $(C_COMPILER_NAME)  $(includes) $(options) $(DIRECTIVES) -c

COMPDEP  = $(CC_COMPILER_NAME) $(includes) $(options) $(cppoptions) $(DIRECTIVES)  -MM -MT
COMPDEPC = $(C_COMPILER_NAME)  $(includes) $(options) $(DIRECTIVES)  -MM -MT

LINK     = $(CC_COMPILER_NAME) $(options) $(linkoptions) $(libs)

######################## SOURCE FILES ##########################

SRC   = $(SRCDIR)inputoutput.cc                         \
        $(SRCDIR)main.cc                                \
        $(SRCDIR)odtline.cc                             \
        $(SRCDIR)processor.cc                           \
        $(SRCDIR)streams.cc                             \
        $(SRCDIR)meshManager.cc                         \
        $(SRCDIR)odtparam.cc                            \
        $(SRCDIR)micromixer.cc                          \
        $(SRCDIR)eddy.cc                                \
        $(SRCDIR)particle.cc                            \
        $(SRCDIR)solver.cc                              \
        $(SRCDIR)kernel.cc                              \
        $(SRCDIR)radiation.cc                           \
        $(SRCDIR)cvodeDriver.cc                         \
        $(SRCDIR)linevariables/lv.cc                    \
        $(SRCDIR)linevariables/lv_pos.cc                \
        $(SRCDIR)linevariables/lv_posf.cc               \
        $(SRCDIR)linevariables/lv_rho.cc                \
        $(SRCDIR)linevariables/lv_rho_const.cc          \
        $(SRCDIR)linevariables/lv_rho_mf.cc             \
        $(SRCDIR)linevariables/lv_dvisc.cc              \
        $(SRCDIR)linevariables/lv_dvisc_const.cc        \
        $(SRCDIR)linevariables/lv_uvw.cc                \
        $(SRCDIR)linevariables/lv_sca.cc                \
        $(SRCDIR)linevariables/lv_enth.cc               \
        $(SRCDIR)linevariables/lv_temp.cc               \
        $(SRCDIR)linevariables/lv_mixf.cc               \
        $(SRCDIR)linevariables/lv_chi.cc                \
        $(SRCDIR)linevariables/lv_hr.cc                 \
        $(SRCDIR)linevariables/lv_ygas.cc               \
        $(SRCDIR)linevariables/lv_ygas_noRxn.cc         \
        $(SRCDIR)linevariables/lv_aDL.cc                \
        $(SRCDIR)linevariables/soot/lv_soot.cc          \
        $(SRCDIR)linevariables/soot/lv_soot_QMOM.cc     \
        $(SRCDIR)linevariables/soot/lv_soot_MONO.cc     \
        $(SRCDIR)linevariables/soot/lv_soot_CQMOM.cc    \
        $(SRCDIR)linevariables/soot/lv_soot_MOMIC.cc    \
        $(SRCDIR)linevariables/soot/CQMOM.cc            \
        $(SRCDIR)linevariables/soot/vandermonde.cc      \
        $(SRCDIR)linevariables/soot/soot_chem.cc        \
        $(SRCDIR)linevariables/soot/pdgen.cc            \
        $(SRCDIR)linevariables/soot/wheeler.cc          \
        $(SRCDIR)linevariables/soot/eispack.cc          \
        $(SRCDIR)nr/interp_linear.cc                    \
        $(SRCDIR)odtcases/odtcase.cc                    \
        $(SRCDIR)odtcases/odtcase_channel.cc            \
        $(SRCDIR)odtcases/odtcase_channelScalar.cc      \
        $(SRCDIR)odtcases/odtcase_RT.cc                 \
        $(SRCDIR)odtcases/odtcase_coldPropaneJet.cc     \
        $(SRCDIR)odtcases/odtcase_coldJet.cc            \
        $(SRCDIR)odtcases/odtcase_jetFlame.cc           \
        $(SRCDIR)odtcases/odtcase_jetMixlRxn.cc         \
        $(SRCDIR)odtcases/odtcase_shearFlow.cc          \
        $(SRCDIR)odtcases/odtcase_HIT.cc        \
        $(SRCDIR)flmlt/odtcase_flmlt.cc         \
        $(SRCDIR)flmlt/lv_ygas_flmlt.cc         \
        $(SRCDIR)flmlt/lv_chi_flmlt.cc          \
        $(SRCDIR)flmlt/lv_enth_flmlt.cc         \
        $(SRCDIR)flmlt/lv_soot_flmlt.cc         \
        $(SRCDIR)flmlt/solver_flmlt.cc          \
        $(SRCDIR)flmlt/micromixer_flmlt.cc      \
        $(SRCDIR)hips/lv_enth_hips.cc           \
        $(SRCDIR)hips/lv_mixf_hips.cc           \
        $(SRCDIR)hips/lv_ygas_hips.cc           \
        $(SRCDIR)hips/lv_ygas_cold_hips.cc      \
        $(SRCDIR)hips/micromixer_hips.cc        \
        $(SRCDIR)hips/solver_hips.cc            \
        $(SRCDIR)hips/odtcase_hips.cc           \
        $(SRCDIR)hips/odtcase_hips_comb.cc      \
        $(SRCDIR)hips/odtcase_hips_simpleRxn.cc 


CVSRC = $(CVODE_SRCDIR)sundials_smalldense.c     \
        $(CVODE_SRCDIR)sundials_nvector.c        \
        $(CVODE_SRCDIR)sundials_math.c           \
        $(CVODE_SRCDIR)sundials_dense.c          \
        $(CVODE_SRCDIR)nvector_serial.c          \
        $(CVODE_SRCDIR)cvode_io.c                \
        $(CVODE_SRCDIR)cvode_dense.c             \
        $(CVODE_SRCDIR)cvode.c

#-------Add your rxn rate file here: make another ifeq block

ifeq ($(CHEMISTRY), ONESTEPC2H4)
    SRC += $(SRCDIR)user_chemical_mechanisms/oneStepRRc2h4.cc
endif
ifeq ($(CHEMISTRY), ONESTEPCH4)
    SRC += $(SRCDIR)user_chemical_mechanisms/oneStepRRch4.cc
endif
ifeq ($(CHEMISTRY), C2H4RED)
    SRC += $(SRCDIR)user_chemical_mechanisms/c2h4RedRR.cc
endif
ifeq ($(CHEMISTRY), CH4RED)
    SRC += $(SRCDIR)user_chemical_mechanisms/ch4RedRR.cc
endif
ifeq ($(CHEMISTRY), SIMPLEDLR)
    SRC += $(SRCDIR)user_chemical_mechanisms/simple_dlr.cc
endif


######################## OBJECT FILES ##########################

OBJS   = $(SRC:.cc=.o)
CVOBJS = $(CVSRC:.c=.o)

######################## BUILD TARGETS ##########################

default : $(PROG)

all : depend default documentation

$(PROG) : $(OBJS) $(CVOBJS)
	$(LINK) $(OBJS) $(CVOBJS) -o $@ $(ext_libs)

%.o: %.cc
	$(COMPCC) $< -o $@
%.o: %.c
	$(COMPC)  $< -o $@

clean :
	rm -f $(OBJS) $(CVOBJS) ./*.mod ./.depend
	rm -f ../run/odt.x
	rm -f ../run/*.log ../run/*.xml
	rm -rf ../doc/doxygen/html/*
	rm -f ../source/user_chemical_mechanisms/*.o

documentation :
	rm -rf ../doc/doxygen/html
	mkdir ../doc/doxygen/html
	cp ../doc/doxygen_build_pages/*.png ../doc/doxygen/html/
	doxygen -s ../doc/Doxyfile




####################### OBJECT DEPENDENCIES #####################

DEPENDS   = $(SRC:.cc=.d)
CVDEPENDS = $(CVSRC:.c=.cvd)

depend :  $(DEPENDS) $(CVDEPENDS)
	cat   $(DEPENDS) $(CVDEPENDS) > .depend
	rm -f $(DEPENDS) $(CVDEPENDS)

%.d :
	$(COMPDEP) $*.o $*.cc > $*.d
%.cvd :
	$(COMPDEPC) $*.o $*.c > $*.cvd

#ifeq ($(wildcard .depend), .depend)
-include .depend
#endif

######################## BUILD YAML  ##########################

yaml :
	./build_yaml.sh

