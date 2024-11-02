SHELL := /bin/bash
PKG := "mbg"
R_EXEC := /usr/bin/R --no-save --quiet

# Build and install package
install:
	@$(R_EXEC) -e "devtools::install()"

# Rebuild man pages
build-docs:
	@$(R_EXEC) -e "devtools::document()"

# Rebuild and deploy man pages
# Requires WEB_USER, WEB_HOST, and WEB_KEY to be exported from parent environment
deploy-docs:
	@$(R_EXEC) -e "devtools::document(); pkgdown::build_site(examples = T, devel = F)"
	@chmod -R o+rX docs/
	@rsync -r -e "ssh -i $(WEB_KEY)" \
        docs/ $(WEB_USER)@$(WEB_HOST):~/public_html/testing/mbg_docs \
       --info=progress2;

# Convenience target to print all of the available targets in this file
# From https://stackoverflow.com/questions/4219255
.PHONY: list
list:
	@LC_ALL=C $(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | \
		awk -v RS= -F: '/^# File/,/^# Finished Make data base/ \
		{if ($$1 !~ "^[#.]") {print $$1}}' | \
		sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'
