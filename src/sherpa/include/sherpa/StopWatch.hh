// 
//  Copyright (C) 2008  Smithsonian Astrophysical Observatory
//
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//

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
