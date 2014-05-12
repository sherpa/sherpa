#ifndef StopWatch_hh
#define StopWatch_hh

#include <ctime>
#include <iostream>

class StopWatch {
public:

  ~StopWatch() {
    std::cout << activity << " took " << get_elapsed_time( ) << " secs\n";
  }

  StopWatch( const std::string& arg ) : activity( arg ), start( std::clock() ){ } 

  double get_elapsed_time( ) const { 
    clock_t total = clock( ) - start; 
    return double( total ) / CLOCKS_PER_SEC ; }

private:

  std::string  activity;
  std::clock_t start;

};

#endif
