default: Acvk

FC = ifort
FFLAGS = 

obj/Q_shep.o: Q_shep.f90
	@mkdir -p obj
	$(FC) $(FFLAGS) -c Q_shep.f90    -o obj/Q_shep.o        -module obj
obj/Acvk.o: Acvk.f90
	$(FC) $(FFLAGS) -c Acvk.f90      -o obj/Acvk.o          -module obj


OBJS =    obj/Q_shep.o          \
          obj/Acvk.o


Acvk: $(OBJS)
	$(FC) $(LNK) -o ~/bin/Acvk $(OBJS) -module obj


clean:
	rm -rf ~/bin/Acvk obj/*.o obj/*.mod

