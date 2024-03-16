.PHONY: formatpy
format-py:
	isort serial
	black serial

.PHONY: formatsh
format-sh:
	find serial -name "*.sh" -exec shellcheck --severity error {} \;
	find serial -name "*.sh" -exec shfmt -w -i 4 {} \;

.PHONY: formatf90
format-f90:
	 find serial -name "*.f90" -exec fprettify {} \;
	
.PHONY: docs
docs:
	doxygen Doxyfile
