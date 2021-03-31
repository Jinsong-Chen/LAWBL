.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}
# .onUnload <- function (libpath) {
#   library.dynam.unload("mypackage", libpath)
# }

StartWelcomeMessage <- function(){

  paste("LAWBL Package ",
        "(version ", utils::packageDescription("LAWBL")$Version,
        "; ",utils::packageDescription("LAWBL")$Date, ")\n",
        "For tutorials, see https://jinsong-chen.github.io/LAWBL/\n",
        sep="")
}
