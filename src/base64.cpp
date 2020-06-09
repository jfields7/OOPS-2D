#include <base64.h>
#include <sstream>

unsigned int Base64::calcSize(unsigned int length){
  unsigned int ret;

  ret = length;
  if(length % 3 != 0){
    ret += 3 - (length % 3);
  }
  ret = (ret/3)*4;

  return ret;
}

std::string Base64::encode(double x){
  static const char charmap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  char *in = reinterpret_cast<char*>(&x);
  unsigned int length = sizeof(double);
  unsigned int end = calcSize(length);

  char *out = new char[end+1];
  out[end] = '\0';
  unsigned int i, j, v;

  for(i = 0, j = 0; i < length; i+=3, j+=4){
    v = in[i];
    v = i+1 < length ? v << 8 | in[i+1] : v << 8;
    v = i+2 < length ? v << 8 | in[i+2] : v << 8;

    out[j] = charmap[(v >> 18) & 0x3F];
    out[j+1] = charmap[(v >> 12) & 0x3F];
    if(i+1 < length){
      out[j+2] = charmap[(v >> 6) & 0x3F];
    }
    else{
      out[j+2] = '=';
    }
    if(i+2 < length){
      out[j+3] = charmap[(v & 0x3F)];
    }
    else{
      out[j+3] = '=';
    }
  }
  std::string str = std::string(out);
  delete[] out;

  return str;
}
