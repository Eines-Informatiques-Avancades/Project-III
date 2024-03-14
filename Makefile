.PHONY: formatpyt
format-py:
	isort src
	black src

.PHONY: formatsh
format-sh:
	find src -name "*.sh" -exec shellcheck --severity error {} \;
	find src -name "*.sh" -exec shfmt -w -i 4 {} \;

.PHONY: formatf90
format-f90:
	 find src -name "*.F90" -exec fprettify {} \;
	
.PHONY: docs
docs:
	doxygen
