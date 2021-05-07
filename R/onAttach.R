.onAttach <- function(libname, pkgname) {
	setup.complete <- (system.file('.perlpath', package = "pipoft") != "")
	if (!setup.complete) {
		warning("Perl setup is incomplete, reconfiguring...")
		prep.setup()
	} 
}