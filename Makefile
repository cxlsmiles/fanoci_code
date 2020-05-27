# FMKMF_F90 not set: using f90
# FMKMF_SFTAG not set: using f90
# FMKMF_SPATH not set: using . 
# FMKMF_LINKOPTS not set: using no link options 
# Main program is $(MAIN).f90 
# process_fsource called with arg $(MAIN).f90 
# $(MAIN).f90 Uses Module globals
# $(MAIN).f90 Uses Module setup_variables
# $(MAIN).f90 Uses Module build_ci_mat
# $(MAIN).f90 Uses Module lapack
# $(MAIN).f90 Uses Module clear_resources
# Full list of modules in $(MAIN).f90: globals setup_variables build_ci_mat lapack clear_resources 
# Uses globals which is in ./globals.f90
# process_fsource called with arg ./globals.f90 
# Full list of modules in ./globals.f90:  
# Uses setup_variables which is in ./setup_variables.f90
# process_fsource called with arg ./setup_variables.f90 
# ./setup_variables.f90 Uses Module globals
# ./setup_variables.f90 Uses Module read_hf_data
# ./setup_variables.f90 Uses Module configurations
# ./setup_variables.f90 Uses Module misc
# ./setup_variables.f90 Uses Module coupling
# Full list of modules in ./setup_variables.f90: globals read_hf_data configurations misc coupling 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses read_hf_data which is in ./read_hf_data.f90
# process_fsource called with arg ./read_hf_data.f90 
# ./read_hf_data.f90 Uses Module globals
# Full list of modules in ./read_hf_data.f90: globals 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses configurations which is in ./configurations.f90
# process_fsource called with arg ./configurations.f90 
# ./configurations.f90 Uses Module globals
# Full list of modules in ./configurations.f90: globals 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses misc which is in ./misc.f90
# process_fsource called with arg ./misc.f90 
# ./misc.f90 Uses Module globals
# Full list of modules in ./misc.f90: globals 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses coupling which is in ./coupling.f90
# process_fsource called with arg ./coupling.f90 
# ./coupling.f90 Uses Module globals
# ./coupling.f90 Uses Module mat_element_el
# ./coupling.f90 Uses Module mat_element_pol
# Full list of modules in ./coupling.f90: globals mat_element_el mat_element_pol 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses mat_element_el which is in ./mat_element_el.f90
# process_fsource called with arg ./mat_element_el.f90 
# ./mat_element_el.f90 Uses Module globals
# ./mat_element_el.f90 Uses Module intindex_module
# Full list of modules in ./mat_element_el.f90: globals intindex_module 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses intindex_module which is in ./intindex_module.f90
# process_fsource called with arg ./intindex_module.f90 
# Full list of modules in ./intindex_module.f90:  
# Uses mat_element_pol which is in ./mat_element_pol.f90
# process_fsource called with arg ./mat_element_pol.f90 
# ./mat_element_pol.f90 Uses Module globals
# ./mat_element_pol.f90 Uses Module intindex_module
# Full list of modules in ./mat_element_pol.f90: globals intindex_module 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses intindex_module which is in ./intindex_module.f90
# ./intindex_module.f90 already in list
# Uses build_ci_mat which is in ./build_ci_mat.f90
# process_fsource called with arg ./build_ci_mat.f90 
# ./build_ci_mat.f90 Uses Module globals
# ./build_ci_mat.f90 Uses Module ci_blocks_el
# ./build_ci_mat.f90 Uses Module ci_blocks_pol
# Full list of modules in ./build_ci_mat.f90: globals ci_blocks_el ci_blocks_pol 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses ci_blocks_el which is in ./ci_blocks_el.f90
# process_fsource called with arg ./ci_blocks_el.f90 
# ./ci_blocks_el.f90 Uses Module globals
# ./ci_blocks_el.f90 Uses Module mat_element_el
# Full list of modules in ./ci_blocks_el.f90: globals mat_element_el 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses mat_element_el which is in ./mat_element_el.f90
# ./mat_element_el.f90 already in list
# Uses ci_blocks_pol which is in ./ci_blocks_pol.f90
# process_fsource called with arg ./ci_blocks_pol.f90 
# ./ci_blocks_pol.f90 Uses Module globals
# ./ci_blocks_pol.f90 Uses Module misc
# ./ci_blocks_pol.f90 Uses Module mat_element_pol
# Full list of modules in ./ci_blocks_pol.f90: globals misc mat_element_pol 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list
# Uses misc which is in ./misc.f90
# ./misc.f90 already in list
# Uses mat_element_pol which is in ./mat_element_pol.f90
# ./mat_element_pol.f90 already in list
# Uses lapack which is in ./lapack.f90
# process_fsource called with arg ./lapack.f90 
# Full list of modules in ./lapack.f90:  
# Uses clear_resources which is in ./clear_resources.f90
# process_fsource called with arg ./clear_resources.f90 
# ./clear_resources.f90 Uses Module globals
# Full list of modules in ./clear_resources.f90: globals 
# Uses globals which is in ./globals.f90
# ./globals.f90 already in list

# ------------------Macro-Defs---------------------
#F90=gfortran
F90=ifort
LIBS = -mkl #/san/local/lib/liblapack++.a 
FLAGS = -i8 -O2
#FLAGS = -O3 -g -traceback -check bounds
# -------------------End-macro-Defs---------------------------

MAIN=fanoci_$(F90)

# Here is the link step 
$(MAIN):globals.o precmod.o stringmod.o read_hf_data.o configurations.o misc.o intindex_module.o reduce_moint.o mat_element_2h.o mat_element_el.o mat_element_pol.o coupling.o setup_variables.o ci_blocks_el.o ci_blocks_pol.o build_ci_mat.o lapack.o clear_resources.o main.o 
	 $(F90) $(FLAGS) -o $(MAIN) globals.o precmod.o stringmod.o read_hf_data.o configurations.o misc.o intindex_module.o reduce_moint.o mat_element_2h.o mat_element_el.o mat_element_pol.o coupling.o setup_variables.o ci_blocks_el.o ci_blocks_pol.o build_ci_mat.o lapack.o clear_resources.o main.o  $(LIBS) 

# Here are the compile steps
 
globals.o:./globals.f90  
	 $(F90) $(FLAGS) -c ./globals.f90 

precmod.o:./precmod.f90
	$(F90) $(FLAGS) -c ./precmod.f90

stringmod.o:./stringmod.f90 precmod.o 
	 $(F90) $(FLAGS) -c ./stringmod.f90

read_hf_data.o:./read_hf_data.f90 globals.o stringmod.o
	 $(F90) $(FLAGS) -c ./read_hf_data.f90 

configurations.o:./configurations.f90 globals.o 
	 $(F90) $(FLAGS) -c ./configurations.f90 

misc.o:./misc.f90 globals.o 
	 $(F90) $(FLAGS) -c ./misc.f90 

intindex_module.o:./intindex_module.f90  
	 $(F90) $(FLAGS) -c ./intindex_module.f90 

reduce_moint.o: ./reduce_moint.f90 globals.o intindex_module.o read_hf_data.o
	$(F90) $(FLAGS) -c ./reduce_moint.f90

mat_element_2h.o:./mat_element_2h.f90 globals.o intindex_module.o 
	 $(F90) $(FLAGS) -c ./mat_element_2h.f90 

mat_element_el.o:./mat_element_el.f90 globals.o intindex_module.o 
	 $(F90) $(FLAGS) -c ./mat_element_el.f90 

mat_element_pol.o:./mat_element_pol.f90 globals.o intindex_module.o 
	 $(F90) $(FLAGS) -c ./mat_element_pol.f90 

coupling.o:./coupling.f90 globals.o mat_element_el.o mat_element_pol.o 
	 $(F90) $(FLAGS) -c ./coupling.f90 

setup_variables.o:./setup_variables.f90 globals.o read_hf_data.o configurations.o misc.o coupling.o intindex_module.o reduce_moint.o
	 $(F90) $(FLAGS) -c ./setup_variables.f90 

ci_blocks_el.o:./ci_blocks_el.f90 globals.o mat_element_el.o 
	 $(F90) $(FLAGS) -c ./ci_blocks_el.f90 

ci_blocks_pol.o:./ci_blocks_pol.f90 globals.o misc.o mat_element_pol.o 
	 $(F90) $(FLAGS) -c ./ci_blocks_pol.f90 

build_ci_mat.o:./build_ci_mat.f90 globals.o ci_blocks_el.o ci_blocks_pol.o mat_element_2h.o
	 $(F90) $(FLAGS) -c ./build_ci_mat.f90 

lapack.o:./lapack.f90  
	 $(F90) $(FLAGS) -c ./lapack.f90 

clear_resources.o:./clear_resources.f90 globals.o 
	 $(F90) $(FLAGS) -c ./clear_resources.f90 

main.o:main.f90 globals.o setup_variables.o build_ci_mat.o lapack.o clear_resources.o 
	 $(F90) $(FLAGS) -c main.f90 
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
  
