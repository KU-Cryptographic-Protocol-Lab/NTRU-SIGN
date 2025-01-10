read("estimate.gp")

cbdMod3(k) = {
  my(p);
  p = floor(2^(2*k)/3);
  [p, p+1, p]/2^(2*k);
};


NTRUPLUS(n,q) = {
  coeffDist = cbdMod3(1);
  Run(n, n, q, coeffDist);
}



run() = {
  print("==================NTRU+SIGN=========================="); 
  print("\n");
  print("=================NTRU+SIGN512========================");  
  NTRUPLUS(512,3329);
  print("\n");
  print("=================NTRU+SIGN1024========================");  
  NTRUPLUS(1024,7681);
  print("\n");
  
}

run();
quit();
