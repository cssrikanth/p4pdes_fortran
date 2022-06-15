include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
FFLAGS += 
LINK.F = $(FLINKER) $(FC_FLAGS) $(FFLAGS) $(FCPPFLAGS) $(LDFLAGS)
e: e.o
	-${LINK.F} -o e e.o  ${PETSC_LIB}
	${RM} e.o

distclean:
	@rm -f *~ e *tmp

