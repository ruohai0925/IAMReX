# AMREX_HOME defines the directory in which we will find AMReX.
AMREX_HOME       ?= ../../../amrex
AMREX_HYDRO_HOME ?= ../../../AMReX-Hydro

TOP = ../..

PRECISION   = DOUBLE
DIM         = 3
COMP        = gnu

DEBUG       = FALSE
USE_MPI     = TRUE
USE_OMP     = FALSE
PROFILE     = FALSE
USE_PARTICLES = TRUE
PARTICLE_PARALLEL = TRUE
USE_CUDA    = FALSE

USE_SENSEI_INSITU = FALSE

EBASE = amr
Blocs := .

include $(TOP)/Exec/Make.IAMR
