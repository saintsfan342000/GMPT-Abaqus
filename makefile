# Names and wildcards
npys := $(wildcard ./ConstructionFiles/*.npy)
#Flags/Paths
.PHONY: clean obj 
clean: obj
	$(RM) $(npys)
obj:
	$(RM) -r  __pycache__

#$(EXEC): $(SRC)
#	$(CC) $(CFLAGS) $(INCFLAGS) $^ $(LDLIBS) $(LDFLAGS) -o $@
