# global variables hack for NSE in data.tables
if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
    'x', 'W',                   # positions
    'xl', 'xr', 'Wl', 'Wr',     # intervals
    'id', 'delta', 'Eu'         # plans
    ))
}

.verbose <- list(
    NONE=0,
    INFO=10,
    TRACE=20,
    DEBUG=30)
