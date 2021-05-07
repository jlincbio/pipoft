require(pipoft)
cat("Hello...testing...\n")

perl.path <- Sys.which("perl")
perl.test <- system.file("p_tests.pl", package = "pipoft")
system2(command = perl.path, args = perl.test)