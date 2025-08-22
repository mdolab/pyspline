%.o: %.F90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

%.o: %.f90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

%.o: %.f
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo

%.o: %.c
	$(CC) $(CC_ALL_FLAGS) -c $< -o $(OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
