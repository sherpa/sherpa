#ifndef STREAMING_EXCEPTION
#define STREAMING_EXCEPTION

//
// Practical C++ Error Handling in Hybrid Environments
// By Gigi Sayfan,    Dr. Dobb's Journal
// Feb 05, 2007
// URL:http://www.ddj.com/cpp/197003350 
//

#include <iostream>
#include <sstream>
#include <memory>
#include <stdexcept>

namespace Sayfan {

  class StreamingException : public std::runtime_error { 
  public:
    StreamingException() : std::runtime_error(""),
			   ss_(std::auto_ptr<std::stringstream>(new std::stringstream())) { }

    ~StreamingException() throw() { }

    template <typename T> StreamingException & operator << (const T & t) {
      (*ss_) << t;
      return *this;
    }
    
    virtual const char * what() const throw() {
      s_ = ss_->str(); return s_.c_str();
    }
    
  private:
    mutable std::auto_ptr<std::stringstream> ss_;
    mutable std::string s_;
  };

  //
  // StreamingException is an exception class that exposes
  // a stream interface and allows one-line formatting of complex error
  // messages. The design and implementation are interesting and use C++
  // constructs such as std::auto_ptr, the mutable modifier, and
  // (begrudgingly) exception specifications.
  //
  // StreamingException is derived from std::runtime_error. It's good form
  // to derive your exception classes from std::standard_error because it
  // communicates your intention (I hope you throw your exceptions at 
  // runtime) and it lets the application catch runtime errors from
  // multiple sources (as long as all the other sources follow this
  // guideline, too). In addition, std::runtime_error provides the
  // what() method that returns a text message with the error description.
  // By default, what() returns the message that was passed in the constructor.
  // StreamingException overrides the virtual what() and returns instead
  // the contents of its stream (it's not called "StreamingException" 
  // for nothing). The StreamingException constructor takes no arguments
  // and passes an empty message to its base std::runtime_error constructor
  // (remember that this message is not used anyway). The constructor then
  // creates a stringstream instance on the heap and puts it in a mutable
  // std::auto_ptr. If you studied your smart pointers well, you know that
  // std::auto_ptr has the unusual property of ownership transfer. When you
  // assign it to another auto_ptr, the assignee gets the ownership of the
  // pointed object and the original object is left with nothing. This peculiar
  // behavior is sometimes dangerous and often gets in the way (you can't put
  // auto_ptr into a standard container because you will lose the ownership to
  // the container). In this case it is exactly what the doctor ordered. You
  // will soon see why and why it must be mutable.
  //
  // The destructor is quite empty, but it can't be dropped. The compiler will
  // indeed generate a default destructor for you, but the default destructor
  // doesn't come with an empty throw() exception specification. This is
  // required because std::runtime_error defines such a virtual destructor.
  // Exception specifications are an annoying misfeature of C++ that specifies
  // what exceptions a method may throw and are part of the method signature.
  // Thankfully, they are optional so you don't see them a lot in the wild. 
  //
  // The templated operator << is where the pedal hits the metal. This is
  // a template method that accepts any streamable T. That translates to
  // standard types plus any other type that implements a conforming
  // operator <<. The implementation simply streams the T argument to its
  // member stringstream s_. The return value is a reference to the
  // StreamingException class itself. This allows chaining multiple calls
  // to operator <<.
  //
  // Finally, StreamingException overrides the virtual what() method to 
  // return the contents of its stringstream. Note the exception
  // specification again. 
  //
  // It's time to reveal StreamingException's secrets. The thrown exception
  // is a temporary object. The instance that is caught in the catch clause is
  // actually a copy of the original exception. I had a couple of choices to
  // keep the error message intact during this copy. Most of them involved
  // implementing a copy constructor and some of them involved duplication
  // of the error message. Instead, I opted to use auto_ptr, which takes care
  // of two issues: It deletes the dynamically allocated stringstream on
  // destruction and it transfers the ownership when StreamingException is
  // copied. Okay, so why mutable? Well, the caught exception is a const
  // reference because the catching code is not supposed to modify the
  // internal state of the exception. However, auto_ptr with its ownership
  // transfer semantics does require a change of state. The mutable modifier
  // was invented exactly for this purpose-being able to modify the internal
  //  state of an object while preserving its conceptual constness.
  //
  
}                                                           // namespace Sayfan

#endif 


/*
#include "StreamingException.hh"

int main (int argc, char * const argv[]) {  
  try {
    if (5 != 3)
      throw StreamingException() << 5 << " is not equal to " << 3;
  } catch (const StreamingException & e) {
    std::cout << e.what() << std::endl;
  }
  return 0;
}
*/
