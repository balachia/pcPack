% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entry-game.R
\name{create.intervals.table}
\alias{create.intervals.table}
\title{Create market information table}
\usage{
create.intervals.table(n)
}
\arguments{
\item{n}{number of intervals to create (n.b. should 1+number of entries)}
}
\description{
market consists of a sequence of enterable intervals
Fields
 id --- (key) entry sequence
 xl, xr --- location of left (right) end; +/-Inf for open intervals
 Wl, Wr --- value at left and right intervals; NA for unentered and open intervals
 state --- interval state
}
