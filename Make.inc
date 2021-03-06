# -*- mode: makefile-gmake -*-
ifneq ($(OS),Windows_NT)
  OS = $(shell uname)
endif

USEGCC = 1
USECLANG = 0

CFLAGS-add += -Wall -lm 

ifeq ($(USECLANG),1)
USEGCC = 0
CC = clang
endif

ifeq ($(USEGCC),1)
CC = g++
endif

ifeq ($(OS), Linux)
SHLIB_EXT = so
CFLAGS-add+=-fPIC
FFLAGS-add+=-fPIC
endif

ifeq ($(OS), FreeBSD)
SHLIB_EXT = so
CFLAGS-add+=-fPIC
FFLAGS-add+=-fPIC
endif

ifeq ($(OS), OpenBSD)
SHLIB_EXT = so
CFLAGS-add+=-fPIC
FFLAGS-add+=-fPIC
endif

ifeq ($(OS), Darwin)
SHLIB_EXT = dylib
CFLAGS-add+=-fPIC
FFLAGS-add+=-fPIC
endif

ifeq ($(OS), Windows_NT)
SHLIB_EXT = dll
endif

# SunOS with gcc
ifeq ($(OS), SunOS)
SHLIB_EXT = so
CFLAGS-add+=-fPIC
FFLAGS-add+=-fPIC
endif