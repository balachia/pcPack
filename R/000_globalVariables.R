# global variables hack for NSE in data.tables
if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
    'x', 'W',                   # positions
    'xl', 'xr', 'Wl', 'Wr',     # intervals
    'id', 'delta', 'Eu'         # plans
    ))
}
