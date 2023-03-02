void reset(){
  cout <<"Resetting ROOT"<<endl;
  gROOT->Reset();
}
void res(){
  reset();
}
void r(){
  reset();
}
string convBase(unsigned long v, long base){
  string digits = "0123456789abcdef";
  string result;
  if((base < 2) || (base > 16)) {
    result = "Error: base out of range.";
  }
  else {
    do {
      result = digits[v % base] + result;
      v /= base;
    }
    while(v);
  }
  return result;
}

