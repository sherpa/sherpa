namespace sherpa {

// 
//  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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


  // Function prototype for 10 parameters.
  template< typename ReType, typename P1=void, typename P2=void, typename P3=void, typename P4=void, typename P5=void, typename P6=void, typename P7=void, typename P8=void, typename P9=void, typename P10=void >
  class FctPtr {
  public:
    typedef ReType (*Fct)( P1, P2, P3, P4, P5, P6, P7, P8, P9, P10 );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator( )( P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6, P7 p7, P8 p8, P9 p9, P10 p10 ) {
      return (*fct)( p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6, typename P7, typename P8, typename P9, typename P10 > inline
  FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10 > fct_ptr( ReType (*fp)( P1, P2, P3, P4, P5, P6, P7, P8, P9, P10 ) ) {
    return FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10 >( fp );
  }

  // Function prototype for 9 parameters.
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6, typename P7, typename P8, typename P9 >
  class FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7, P8, P9, void > {
  public:
    typedef ReType (*Fct)( P1, P2, P3, P4, P5, P6, P7, P8, P9 );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator( )( P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6, P7 p7, P8 p8, P9 p9 ) {
      return (*fct)( p1, p2, p3, p4, p5, p6, p7, p8, p9 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6, typename P7, typename P8, typename P9 > inline
  FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7, P8, P9 > fct_ptr( ReType (*fp)( P1, P2, P3, P4, P5, P6, P7, P8, P9 ) ) {
    return FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7, P8, P9 >( fp );
  }

  // Function prototype for 8 parameters.
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6, typename P7, typename P8 >
  class FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7, P8, void > {
  public:
    typedef ReType (*Fct)( P1, P2, P3, P4, P5, P6, P7, P8 );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator( )( P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6, P7 p7, P8 p8 ) {
      return (*fct)( p1, p2, p3, p4, p5, p6, p7, p8 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6, typename P7, typename P8 > inline
  FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7, P8 > fct_ptr( ReType (*fp)( P1, P2, P3, P4, P5, P6, P7, P8 ) ) {
    return FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7, P8 >( fp );
  }

  // Function prototype for 7 parameters.
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6, typename P7 >
  class FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7, void > {
  public:
    typedef ReType (*Fct)( P1, P2, P3, P4, P5, P6, P7 );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator( )( P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6, P7 p7 ) {
      return (*fct)( p1, p2, p3, p4, p5, p6, p7 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6, typename P7 > inline
  FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7 > fct_ptr( ReType (*fp)( P1, P2, P3, P4, P5, P6, P7 ) ) {
    return FctPtr< ReType, P1, P2, P3, P4, P5, P6, P7 >( fp );
  }

  // Function prototype for 6 parameters.
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6 >
  class FctPtr< ReType, P1, P2, P3, P4, P5, P6, void > {
  public:
    typedef ReType (*Fct)( P1, P2, P3, P4, P5, P6 );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator( )( P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6 ) {
      return (*fct)( p1, p2, p3, p4, p5, p6 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6 > inline
  FctPtr< ReType, P1, P2, P3, P4, P5, P6 > fct_ptr( ReType (*fp)( P1, P2, P3, P4, P5, P6 ) ) {
    return FctPtr< ReType, P1, P2, P3, P4, P5, P6 >( fp );
  }

  // Function prototype for 5 parameters.
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5 >
  class FctPtr< ReType, P1, P2, P3, P4, P5, void > {
  public:
    typedef ReType (*Fct)( P1, P2, P3, P4, P5 );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator( )( P1 p1, P2 p2, P3 p3, P4 p4, P5 p5 ) {
      return (*fct)( p1, p2, p3, p4, p5 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1, typename P2, typename P3, typename P4, typename P5 > inline
  FctPtr< ReType, P1, P2, P3, P4, P5 > fct_ptr( ReType (*fp)( P1, P2, P3, P4, P5 ) ) {
    return FctPtr< ReType, P1, P2, P3, P4, P5 >( fp );
  }

  // Function prototype for 4 parameters.
  template<typename ReType, typename P1, typename P2, typename P3, typename P4>
  class FctPtr< ReType, P1, P2, P3, P4, void > {
  public:
    typedef ReType (*Fct)( P1, P2, P3, P4 );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator( )( P1 p1, P2 p2, P3 p3, P4 p4 ) {
      return (*fct)( p1, p2, p3, p4 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1, typename P2, typename P3, typename P4 > inline
  FctPtr< ReType, P1, P2, P3, P4 > fct_ptr( ReType (*fp)( P1, P2, P3, P4 ) ) {
    return FctPtr< ReType, P1, P2, P3, P4 >( fp );
  }

  // Function prototype for 3 parameters.
  template< typename ReType, typename P1, typename P2, typename P3 >
  class FctPtr< ReType, P1, P2, P3, void, void > {
  public:
    typedef ReType (*Fct)( P1, P2, P3 );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator( )( P1 p1, P2 p2, P3 p3 ) { 
      return (*fct)( p1, p2, p3 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1, typename P2, typename P3 > inline
  FctPtr< ReType, P1, P2, P3 > fct_ptr( ReType (*fp)( P1, P2, P3 ) ) {
    return FctPtr< ReType, P1, P2, P3 >( fp );
  }

  // Function prototype for 2 parameters.
  template< typename ReType, typename P1, typename P2 >
  class FctPtr< ReType, P1, P2, void, void, void > {
  public:
    typedef ReType (*Fct)( P1, P2 );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator( )( P1 p1, P2 p2 ) { return (*fct)( p1, p2 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1, typename P2 > inline
  FctPtr< ReType, P1, P2 > fct_ptr( ReType (*fp)( P1, P2 ) ) {
    return FctPtr< ReType, P1, P2 >( fp );
  }

  // Function prototype for 1 parameter.
  template< typename ReType, typename P1 >
  class FctPtr< ReType, P1, void, void, void, void > {
  public:
    typedef ReType (*Fct)( P1 );
    FctPtr( Fct f, P1 p1 ) : fct( f ) { }
    ReType operator() ( P1 p1 ) { return (*fct)( p1 ); }
  private:
    Fct fct;
  };
  template< typename ReType, typename P1 > inline
  FctPtr< ReType, P1 > fct_ptr( ReType (*fp)( P1 ) ) {
    return FctPtr< ReType, P1 >( fp );
  }

  // Function prototype for no parameter.
  template< typename ReType >
  class FctPtr< ReType, void, void, void, void, void > {
  public:
    typedef ReType (*Fct)( );
    FctPtr( Fct f ) : fct( f ) { }
    ReType operator() ( ) { return (*fct)( ); }
  private:
    Fct fct;
  };
  template< typename ReType > inline
  FctPtr< ReType > fct_ptr( ReType (*fp)( ) ) {
    return FctPtr< ReType >( fp );
  }

}                                                           // namespace sherpa

/*
#include <iostream>

int hi0( ) {
  std::cout << "hi0()\n";
  return 0;
}

int hi1( int& arg1 ) {
  std::cout << "hi1(" << arg1 << ")\n";
  arg1 = - arg1;
  return arg1;
}

int hi2( int arg1, int arg2 ) {
  std::cout << "hi2(" << arg1 << "," << arg2 << ")\n";
  return arg1 + arg2;
}

void hi22( int arg1, int* arg2 ) {
  std::cout << "hi22( " << arg1 << ", " << *arg2 << ")\n";
}

int hi3( int arg1, int arg2, int arg3 ) {
  std::cout << "hi3(" << arg1 << "," << arg2 << "," << arg3 << ")\n";
  return arg1 + arg2 + arg3;
}

int hi4( int arg1, int arg2, int arg3, int arg4 ) {
  std::cout << "hi4(" << arg1 << "," << arg2 << "," << arg3 << ',' << arg4 << ")\n";
  return arg1 + arg2 + arg3 + arg4;
}

int hi5(int arg1, int arg2, int arg3, int arg4, int arg5 ) {
  std::cout << "hi5(" << arg1 << "," << arg2 << "," << arg3 << "," << arg4 << "," << arg5 << ")\n";
  return arg1 + arg2 + arg3 + arg4 + arg5;
}

int hi6(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6 ) {
  std::cout << "hi6(" << arg1 << "," << arg2 << "," << arg3 << "," << arg4 << "," << arg5 << "," << arg6 << ")\n";
  return arg1 + arg2 + arg3 + arg4 + arg5 + arg6;
}

int hi7(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7 ) {
  std::cout << "hi7(" << arg1 << "," << arg2 << "," << arg3 << "," << arg4 << "," << arg5 << "," << arg6 << "," << arg7 << ")\n";
  return arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7;
}

int hi8(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8 ) {
  std::cout << "hi8(" << arg1 << "," << arg2 << "," << arg3 << "," << arg4 << "," << arg5 << "," << arg6 << "," << arg7 << "," << arg8 << ")\n";
  return arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8;
}

int hi9(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9 ) {
  std::cout << "hi9(" << arg1 << "," << arg2 << "," << arg3 << "," << arg4 << "," << arg5 << "," << arg6 << "," << arg7 << "," << arg8 << "," << arg9 << ")\n";
  return arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8 + arg9;
}

int hi10(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9, int arg10 ) {
  std::cout << "hi10(" << arg1 << "," << arg2 << "," << arg3 << "," << arg4 << "," << arg5 << "," << arg6 << "," << arg7 << "," << arg8 << "," << arg9 << "," << arg10 << ")\n";
  return arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8 + arg9 + arg10;
}

template < typename Functor >
void demo( Functor func ) {

  std::cout << "demo() : " << func() << '\n';

}

int main( int argc, char* arg[] ) {

  int arg1 = 1, arg2=2, arg3=3, arg4=4, arg5=5, arg6=6, arg7=7, arg8=8, arg9=9,
    arg10=10;

  {
    sherpa::FctPtr<int> func0( hi0 );
    func0( );
    // note that using fct_ptr then one does not have to declare prototype!
    demo( sherpa::fct_ptr( hi0 ) );
  }

  sherpa::FctPtr<int,int&> func1( hi1, arg1 );
  func1( arg1 );
  std::cout << "arg1 = " << arg1 << '\n';

  sherpa::FctPtr<int,int,int> func2( hi2 );
  func2( arg1, arg2 );

  sherpa::FctPtr<void,int,int*> func22( hi22 );
  func22( arg1, &arg2 );

  sherpa::FctPtr<int,int,int,int> func3( hi3 );
  func3( arg1, arg2, arg3 );

  sherpa::FctPtr<int,int,int,int,int> func4( hi4 );
  func4( arg1, arg2, arg3, arg4 );

  sherpa::FctPtr<int,int,int,int,int,int> func5( hi5 );
  func5( arg1, arg2, arg3, arg4, arg5 );

  sherpa::FctPtr<int,int,int,int,int,int,int> func6( hi6 );
  func6( arg1, arg2, arg3, arg4, arg5, arg6 );

  sherpa::FctPtr<int,int,int,int,int,int,int,int> func7( hi7 );
  func7( arg1, arg2, arg3, arg4, arg5, arg6, arg7 );

  sherpa::FctPtr<int,int,int,int,int,int,int,int,int> func8( hi8 );
  func8( arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8 );

  sherpa::FctPtr<int,int,int,int,int,int,int,int,int,int> func9( hi9 );
  func9( arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9 );

  sherpa::FctPtr<int,int,int,int,int,int,int,int,int,int,int> func10( hi10 );
  func10( arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10 );

  return 0;

}
*/
