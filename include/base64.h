#ifndef BASE64_H
#define BASE64_H
#include <string>

/**
 * A class to help with encoding data in Base64, primarily for VTK binary output.
 * Based on the code found at 
 * https://nachtimwald.com/2017/11/18/base64-encode-and-decode-in-c/
 */
class Base64{
  private:
    static unsigned int calcSize(unsigned int length);
  public:
    static std::string encode(double x);
};

#endif
