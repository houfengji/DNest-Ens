#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>

class Exception {   // How to throw an exception
private: 
	std::string MessageText;
public:
	Exception(const std::string & ErrorMessage);  //  Set the exception message 
	std::string ExceptionMessage();               //  Retrieve the exception message outside the construcdtor
};

#endif
