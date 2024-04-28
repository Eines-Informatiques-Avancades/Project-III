.PHONY: formatf90
format-f90:
	find parallel -name "*.f90" -exec fprettify {} \;

.PHONY: docs
docs:
	doxygen Doxyfile
	cd docs/latex && make pdf && cd ..
