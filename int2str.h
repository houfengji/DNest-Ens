#ifndef INT2STR_H
#define INT2STR_H

#include <string>
void int2str( const unsigned long i,
              const size_t width,   // how many digits to output
              std::string & str);

std::string int2str(const unsigned long i,
                    const size_t width);   // how many digits to output

std::string int2str(const unsigned long i); 

#endif
