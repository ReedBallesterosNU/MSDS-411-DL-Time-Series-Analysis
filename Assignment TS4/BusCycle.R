
# By Mark Heiple

# This function implements the k = … equation from Tsay p. 42 (the one with a and b, not phi).
#c is a complex number

period <- function(c) {

  #get real and imaginary components
  r = Re(c)
  j = Im(c)

  #magnitude
  m = sqrt(r^2 + j^2)
  period = (2 * pi)/acos(r/m)

  return( period )
}

#extract the complex roots from list returned from polyroot()

all_complex = function(x) {
  i = round(Im(x),digits=3) != 0
  return(x[i])
}
