#include <vector>

/*
 *
 * Define a map (lookup table) that contains boundaries for binned data
 * in two dimensions. 
 *
 * Functions allow to return the bin # for an arbitrary point within the 
 * allowed phase space.
 *
 */

class BinnedData2D
{
public:
  typedef std::pair<double,double> dpair_t;
  
  class Boundary2D : public std::pair<dpair_t, dpair_t>
  {
  public:
    Boundary2D( double xmin, double xmax, double ymin, double ymax ) : std::pair<dpair_t, dpair_t>( dpair_t( xmin, xmax ),
												    dpair_t( ymin, ymax ) ) { }
    
    double xmin() const { return this->first.first; }
    double xmax() const { return this->first.second; }
    double ymin() const { return this->second.first; }
    double ymax() const { return this->second.second; }

    bool contains( double x, double y ) const { return ( x <= xmax() && x >= xmin() && y <= ymax() && y >= ymin() ); }
    bool overlaps( const Boundary2D& bndry ) const;
  };

public:

  BinnedData2D( double xmin, double xmax, double ymin, double ymax );
  ~BinnedData2D() {}

  // add a bin to the lookup table
  //    xmin, xmax: boundaries in x dimension
  //    ymin, ymax: boundaries in y dimension
  //
  // returns an integer id for this bin, -1 if failed to add bin
  int add_bin( double xmin, double xmax, double ymin, double ymax );
  
  // find bin that contains this point
  //    x: location along x axis
  //    y: location along y axis
  //
  // returns the integer id for bin containing this point, -1 if failed to find bin
  int get_bin( double x, double y );

private:
  
  std::vector< Boundary2D > _bins;

  double _xmin;
  double _xmax; 
  double _ymin;
  double _ymax;
};
  
BinnedData2D::BinnedData2D( double xmin, double xmax, double ymin, double ymax ) 
  : _bins( ) 
{
  _xmin = xmin;
  _xmax = xmax;
  _ymin = ymin;
  _ymax = ymax;
}

int BinnedData2D::add_bin( double xmin, double xmax, double ymin, double ymax )
{
  if( xmin < _xmin || xmax > _xmax || ymin < _ymin || ymax > _ymax )
    return -1;

  if( xmin >= xmax || ymin >= ymax )
    return -1;

  Boundary2D bin( xmin, xmax, ymin, ymax );
  
  for( unsigned int ibin=0; ibin < _bins.size(); ++ibin )
  {
    if( _bins[ibin].overlaps( bin ) )
      return -1;
  }
  
  _bins.push_back( bin );
  return _bins.size();
}

int BinnedData2D::get_bin( double x, double y )
{
  for( unsigned int ibin=0; ibin < _bins.size(); ++ibin )
  {
    if( _bins[ibin].contains( x, y ) )
      return ibin;
  }
  return -1;
}

bool BinnedData2D::Boundary2D::overlaps( const Boundary2D& bndry ) const
{
  return !( bndry.xmin() >= xmax() ||
	    bndry.xmax() <= xmin() || 
	    bndry.ymax() <= ymin() ||
	    bndry.ymin() >= ymax() );
}
