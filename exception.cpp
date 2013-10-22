#include <string>
#include "exception.h"

Exception::Exception(const std::string & ErrorMessage){  //  Set the exception message
	MessageText = ErrorMessage;
}

std::string Exception::ExceptionMessage(){
	return MessageText;
}
