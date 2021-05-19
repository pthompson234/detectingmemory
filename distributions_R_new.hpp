template<class Type>
Type dvonmises(Type y, Type mu, Type kappa, int give_log = 0) {
  Type ret = 0;
  
  if (kappa == 0) {
    ret = log(1/(2 * M_PI));
  } else if (kappa < 100000) {
    ret = kappa * (cos(y - mu)) - (log(2 * M_PI) + log(besselI(kappa, Type(0))));
  } else if (y == mu) {
    ret = Type(-INFINITY); // not exact but couldn't figure the modulus operator out here
  }
  
  return ( give_log ? ret : exp(ret) );
}
VECTORIZE4_ttti(dvonmises)
