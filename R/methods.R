TransMethod <- function(object, s, t, ...) {
	UseMethod("TransMethod");
} # TransMethod

TransBoot <- function(object, s, t, ...) {
	UseMethod("TransBoot");
} # TransBoot

TPBoot <- function(object, UT, nboot, ...) {
	UseMethod("TPBoot");
} # TPBoot

TPCBoot <- function(object, UT, UX, nboot, ...) {
	UseMethod("TPCBoot");
} # TPCBoot

TransPROB <- function(object, UT, ...) {
	UseMethod("TransPROB");
} # TransPROB

toTPmsm <- function(lst, UT, s, t, statenames) {
	UseMethod("toTPmsm");
} # toTPmsm

BtoTPmsm <- function(lst, UT, s, t, statenames, nboot, conflevel, methodboot) {
	UseMethod("BtoTPmsm");
} # BtoTPmsm

TransMatrix <- function(x) {
	UseMethod("TransMatrix");
} # TransMatrix
