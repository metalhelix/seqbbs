

BUILT=$(shell ls -t *.tar.gz | head -1)

all: test

test: build
	@echo MY_VAR IS $(BUILT)
	R CMD check $(BUILT)

build:
	R CMD build .
